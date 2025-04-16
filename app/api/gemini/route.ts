import { NextRequest } from 'next/server';
import { loadContext, getSystemPromptWithContext } from '@/app/utils/context';
import { getProgramById, loadConfig } from '@/app/utils/config';
import { generateGeminiStreamingResponse } from '@/app/utils/gemini-client';

// This is a streaming version of the Gemini API
export const runtime = "edge";

export async function POST(request: NextRequest) {
  try {
    // Check if API key is configured
    if (!process.env.GEMINI_API_KEY) {
      return new Response(
        JSON.stringify({
          error: 'Gemini API key is not configured',
          details: 'Please set the GEMINI_API_KEY environment variable'
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

    // Load context for the program
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

    // Log API key for debugging (partial key only)
    const apiKeyPrefix = process.env.GEMINI_API_KEY?.substring(0, 10);
    console.log(`Using Gemini API key starting with: ${apiKeyPrefix}...`);

    // Create system message with context
    const systemPrompt = getSystemPromptWithContext(programId, context);

    // Import the parseModelId function
    const { parseModelId } = await import('@/app/utils/config');

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

    // Log the model selection
    console.log('Using Gemini model:', modelName, modelId ? `(from model ID: ${modelId})` : '(from default config)');

    // Log the request for debugging
    console.log('Gemini API Request:', {
      programId,
      model: modelName,
      messageCount: messages.length,
      systemPromptLength: systemPrompt.length
    });

    try {
      // Log the model options being used
      console.log('Using model options:', modelOptions);

      // Generate a streaming response
      const stream = await generateGeminiStreamingResponse(
        process.env.GEMINI_API_KEY,
        modelName,
        systemPrompt,
        messages.map((msg: any) => ({
          role: msg.role,
          content: msg.content
        })),
        modelOptions.temperature ?? 0.7, // Use model-specific temperature or default
        modelOptions.max_completion_tokens // Use model-specific max tokens if available
      );

      console.log('Gemini stream created successfully');

      // Create a response with the stream
      return new Response(stream, {
        headers: {
          'Content-Type': 'text/event-stream',
          'Cache-Control': 'no-cache',
          'Connection': 'keep-alive'
        }
      });
    } catch (error) {
      // Handle model-specific errors
      console.error('Error creating Gemini completion:', error);

      // Log detailed error information
      const errorMessage = error instanceof Error ? error.message : String(error);
      const errorName = error instanceof Error ? error.name : 'Unknown';
      console.error('Detailed error info:', {
        name: errorName,
        message: errorMessage,
        stack: error instanceof Error ? error.stack : 'No stack trace',
        model: modelName,
        apiKeyConfigured: !!process.env.GEMINI_API_KEY
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
            error: 'Gemini API key error',
            details: 'Please check your Gemini API key configuration.',
            message: errorMessage
          }),
          { status: 401, headers: { 'Content-Type': 'application/json' } }
        );
      }

      // Re-throw for general error handling
      throw error;
    }
  } catch (error) {
    console.error('Error in Gemini API:', error);

    // Create a more detailed error message
    const errorMessage = error instanceof Error
      ? `Error: ${error.message}`
      : 'An unknown error occurred';

    // Log additional context for debugging
    console.error('Error context:', {
      programId: request.body ? 'Present' : 'Missing',
      geminiApiKey: process.env.GEMINI_API_KEY ? 'Configured' : 'Missing',
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
