import { NextRequest } from 'next/server';
import OpenAI from 'openai';
import { loadContext, getSystemPromptWithContext } from '@/app/utils/context';
import { getProgramById } from '@/app/utils/config';

// This enables Edge Functions in Vercel
export const runtime = "edge";

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
      apiKey: process.env.OPENAI_API_KEY,
      dangerouslyAllowBrowser: true // Required for edge runtime
    });

    // Log API key for debugging (partial key only)
    const apiKeyPrefix = process.env.OPENAI_API_KEY?.substring(0, 10);
    console.log(`Using API key starting with: ${apiKeyPrefix}...`);

    // Create system message with context
    const systemMessage = {
      role: "system",
      content: getSystemPromptWithContext(programId, context)
    };

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
    console.log('OpenAI API Request:', {
      programId,
      model: modelName,
      messageCount: messages.length,
      systemMessageLength: systemMessage.content.length
    });

    // Log the request for debugging
    console.log('Using OpenAI model:', modelName, modelId ? `(from model ID: ${modelId})` : '(from default config)');

    try {
      // Log the model options being used
      console.log('Using model options:', modelOptions);

      // Create the completion options
      const completionOptions: any = {
        model: modelName,
        messages: [
          systemMessage,
          ...messages.map((msg: any) => ({
            role: msg.role,
            content: msg.content
          }))
        ],
        stream: true
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

      // Add any other model-specific options
      Object.entries(modelOptions)
        .filter(([key]) => !['temperature', 'max_completion_tokens'].includes(key))
        .forEach(([key, value]) => {
          completionOptions[key] = value;
        });

      // Make sure stream is set to true
      completionOptions.stream = true;

      // Create the stream - use the correct type for streaming responses
      const stream = await openai.chat.completions.create(completionOptions) as any;

      console.log('Stream created successfully');

      // Create a response with the stream
      const response = new Response(stream);

      // Log response headers for debugging
      console.log('Response created with headers:', {
        'content-type': response.headers.get('content-type'),
        'transfer-encoding': response.headers.get('transfer-encoding')
      });

      return response;
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

    // This line is unreachable due to the try/catch above
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
