import { NextRequest } from 'next/server';
import OpenAI from 'openai';
import { loadContext, getSystemPromptWithContext } from '@/app/utils/context';
import { getProgramById } from '@/app/utils/config';

// This is a non-streaming version of the chat API
// It's used for better compatibility and simpler implementation

export async function POST(request: NextRequest) {
  try {
    // Check if API key is configured
    if (!process.env.OPENAI_API_KEY) {
      return new Response(
        JSON.stringify({
          error: 'OpenAI API key is not configured',
          details: 'Please set the OPENAI_API_KEY environment variable'
        }),
        { status: 500, headers: { 'Content-Type': 'application/json' } }
      );
    }

    // Parse message from post
    const { programId, messages, modelId } = await request.json();

    // Get program configuration
    const program = getProgramById(programId);
    if (!program) {
      return new Response(
        JSON.stringify({ error: `Program ${programId} not found` }),
        { status: 400 }
      );
    }

    // Load context for the program (now returns a Promise)
    let contextToLoad: string | string[] = [];

    // Check if program has contextFiles array
    if (program && 'contextFiles' in program && Array.isArray(program.contextFiles)) {
      contextToLoad = program.contextFiles;
    }
    // Fallback to single contextFile if available
    else if (program && 'contextFile' in program && typeof program.contextFile === 'string') {
      contextToLoad = program.contextFile;
    }

    const context = await loadContext(contextToLoad);

    // Create OpenAI client with explicit configuration
    const openai = new OpenAI({
      apiKey: process.env.OPENAI_API_KEY
    });

    // Log API key for debugging (partial key only)
    const apiKeyPrefix = process.env.OPENAI_API_KEY?.substring(0, 10);
    console.log(`Using API key starting with: ${apiKeyPrefix}...`);

    // Create system message with context
    const systemMessage = {
      role: "system",
      content: getSystemPromptWithContext(programId, context)
    };

    // Log the request for debugging
    console.log('OpenAI API Request:', {
      programId,
      model: process.env.OPENAI_MODEL || "gpt-4.1-mini-2025-04-14",
      messageCount: messages.length,
      systemMessageLength: systemMessage.content.length
    });

    // Import the parseModelId and loadConfig functions
    const { parseModelId, loadConfig } = await import('@/app/utils/config');

    // Get the model name from the model ID or use default from config
    const config = loadConfig();
    let modelName: string;
    let modelOptions: any = {};

    // If a model ID is provided, parse it to get the model name and options
    if (modelId) {
      try {
        const { modelName: parsedModelName } = parseModelId(modelId);
        modelName = parsedModelName;

        // Find the model in the config to get its options
        const modelConfig = config.availableModels.find(model => model.id === modelId);
        if (modelConfig?.options) {
          modelOptions = modelConfig.options;
        }
      } catch (error) {
        console.warn(`Invalid model ID format: ${modelId}. Using default model.`);
        // Use default model from config
        const { modelName: defaultModelName } = parseModelId(config.defaultModelId);
        modelName = defaultModelName;

        // Get options for the default model
        const defaultModelConfig = config.availableModels.find(model => model.id === config.defaultModelId);
        if (defaultModelConfig?.options) {
          modelOptions = defaultModelConfig.options;
        }
      }
    } else {
      // Use default model from config
      const { modelName: defaultModelName } = parseModelId(config.defaultModelId);
      modelName = defaultModelName;

      // Get options for the default model
      const defaultModelConfig = config.availableModels.find(model => model.id === config.defaultModelId);
      if (defaultModelConfig?.options) {
        modelOptions = defaultModelConfig.options;
      }
    }

    // Log the request for debugging
    console.log('Using OpenAI model:', modelName, modelId ? `(from model ID: ${modelId})` : '(from default config)');

    // Log the model options being used
    console.log('Using model options:', modelOptions);

    try {
      // Create the completion options
      const completionOptions: any = {
        model: modelName,
        messages: [
          systemMessage,
          ...messages.map((msg: any) => ({
            role: msg.role,
            content: msg.content
          }))
        ]
      };

      // Add temperature if supported by the model
      if (modelOptions.temperature !== undefined) {
        completionOptions.temperature = modelOptions.temperature;
      }

      // Add max tokens if specified - handle different parameter names for different models
      if (modelOptions.max_completion_tokens) {
        // o4-mini uses max_completion_tokens, other models use max_tokens
        if (modelName === 'o4-mini') {
          completionOptions.max_completion_tokens = modelOptions.max_completion_tokens;
        } else {
          completionOptions.max_tokens = modelOptions.max_completion_tokens;
        }
      }

      // Create a non-streaming completion
      const completion = await openai.chat.completions.create(completionOptions);

      // Log token usage
      console.log('Completion created successfully');
      console.log('Token usage:', {
        promptTokens: completion.usage?.prompt_tokens || 0,
        completionTokens: completion.usage?.completion_tokens || 0,
        totalTokens: completion.usage?.total_tokens || 0,
        model: modelName,
        programId
      });

      // Return the completion
      return new Response(
        JSON.stringify({
          content: completion.choices[0].message.content,
          role: 'assistant'
        }),
        {
          status: 200,
          headers: { 'Content-Type': 'application/json' }
        }
      );
    } catch (error) {
      // Handle model-specific errors
      console.error('Error creating chat completion:', error);

      // Log detailed error information
      const errorMessage = error instanceof Error ? error.message : String(error);
      const errorName = error instanceof Error ? error.name : 'Unknown';
      console.error('Detailed error info:', {
        name: errorName,
        message: errorMessage,
        stack: error instanceof Error ? error.stack : 'No stack trace',
        model: modelName,
        apiKeyConfigured: !!process.env.OPENAI_API_KEY
      });

      // Check if it's a model-related error
      if (errorMessage.includes('model')) {
        return new Response(
          JSON.stringify({
            error: 'Invalid model configuration',
            details: `The specified model "${modelName}" may not be available or properly configured.`,
            message: errorMessage
          }),
          { status: 400, headers: { 'Content-Type': 'application/json' } }
        );
      }

      // Check if it's an API key error
      if (errorMessage.includes('API key') || errorMessage.includes('authentication')) {
        return new Response(
          JSON.stringify({
            error: 'OpenAI API key error',
            details: 'Please check your OpenAI API key configuration.',
            message: errorMessage
          }),
          { status: 401, headers: { 'Content-Type': 'application/json' } }
        );
      }

      // Re-throw for general error handling
      throw error;
    }
  } catch (error) {
    console.error('Error in chat API:', error);

    // Create a more detailed error message
    const errorMessage = error instanceof Error
      ? `Error: ${error.message}`
      : 'An unknown error occurred';

    // Log additional context for debugging
    console.error('Error context:', {
      programId: request.body ? 'Present' : 'Missing',
      openaiApiKey: process.env.OPENAI_API_KEY ? 'Configured' : 'Missing',
      errorType: error instanceof Error ? error.constructor.name : typeof error,
      timestamp: new Date().toISOString()
    });

    // Return a more informative error response
    return new Response(
      JSON.stringify({
        error: 'An error occurred processing your request',
        details: errorMessage,
        timestamp: new Date().toISOString()
      }),
      {
        status: 500,
        headers: {
          'Content-Type': 'application/json'
        }
      }
    );
  }
}
