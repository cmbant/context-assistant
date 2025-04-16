import { GoogleGenAI } from '@google/genai';

// Cache for storing AI client instances
const clientCache: Record<string, GoogleGenAI> = {};

/**
 * Get a Google Gen AI client instance
 * @param apiKey The Google API key
 * @returns A Google Gen AI client instance
 */
export function getGoogleGenAI(apiKey: string): GoogleGenAI {
  // Check if the client is already cached
  if (clientCache[apiKey]) {
    return clientCache[apiKey];
  }

  // Create a new Google Gen AI instance
  const ai = new GoogleGenAI({ apiKey });

  // Cache the client
  clientCache[apiKey] = ai;

  return ai;
}

/**
 * Generate a response from Gemini
 * @param apiKey The Google API key
 * @param modelName The model name (e.g., 'gemini-2.0-flash')
 * @param systemPrompt The system prompt
 * @param messages The conversation history
 * @param temperature The temperature for generation
 * @param programId Optional program ID for token usage logging
 * @param maxOutputTokens Optional maximum number of tokens to generate
 * @returns The generated response
 */
export async function generateGeminiResponse(
  apiKey: string,
  modelName: string,
  systemPrompt: string,
  messages: { role: string; content: string }[],
  temperature: number = 0.7,
  programId?: string,
  maxOutputTokens?: number
): Promise<string> {
  try {
    // Get the AI client
    const ai = getGoogleGenAI(apiKey);

    // Create a chat session
    const chat = ai.chats.create({
      model: modelName,
      history: [
        {
          role: "user",
          parts: [{ text: systemPrompt }]
        },
        {
          role: "model",
          parts: [{ text: "I understand. I will follow these instructions." }]
        }
      ]
    });

    // Add the conversation history
    for (let i = 0; i < messages.length - 1; i++) {
      const message = messages[i];
      // Convert our message format to Gemini's format
      if (message.role === 'assistant') {
        // For assistant messages, add as model response
        await chat.sendMessage({
          message: message.content
        });
      } else {
        // For user messages
        await chat.sendMessage({
          message: message.content
        });
      }
    }

    // Send the last message and get the response
    const response = await chat.sendMessage({
      message: messages[messages.length - 1].content,
      config: {
        temperature,
        maxOutputTokens: maxOutputTokens || 4096, // Use provided max tokens or default
      }
    });

    // Estimate token usage (Gemini doesn't provide token counts directly)
    const promptText = [
      systemPrompt,
      ...messages.map(m => m.content)
    ].join(' ');
    const responseText = response.text || '';

    // Rough estimation: ~4 chars per token for English text
    const estimatedPromptTokens = Math.ceil(promptText.length / 4);
    const estimatedResponseTokens = Math.ceil(responseText.length / 4);

    // Log estimated token usage
    console.log('Estimated Gemini token usage:', {
      promptTokens: estimatedPromptTokens,
      completionTokens: estimatedResponseTokens,
      totalTokens: estimatedPromptTokens + estimatedResponseTokens,
      model: modelName,
      programId: programId || 'unknown'
    });

    return response.text || '';
  } catch (error) {
    console.error('Error generating Gemini response:', error);
    throw error;
  }
}

/**
 * Generate a streaming response from Gemini
 * @param apiKey The Google API key
 * @param modelName The model name (e.g., 'gemini-2.0-flash')
 * @param systemPrompt The system prompt
 * @param messages The conversation history
 * @param temperature The temperature for generation
 * @param maxOutputTokens Optional maximum number of tokens to generate
 * @returns A ReadableStream of the generated response
 */
export async function generateGeminiStreamingResponse(
  apiKey: string,
  modelName: string,
  systemPrompt: string,
  messages: { role: string; content: string }[],
  temperature: number = 0.7,
  maxOutputTokens?: number
): Promise<ReadableStream> {
  try {
    // Get the AI client
    const ai = getGoogleGenAI(apiKey);

    // Create a chat session
    const chat = ai.chats.create({
      model: modelName,
      history: [
        {
          role: "user",
          parts: [{ text: systemPrompt }]
        },
        {
          role: "model",
          parts: [{ text: "I understand. I will follow these instructions." }]
        }
      ]
    });

    // Add the conversation history
    for (let i = 0; i < messages.length - 1; i++) {
      const message = messages[i];
      // Convert our message format to Gemini's format
      if (message.role === 'assistant') {
        // For assistant messages, add as model response
        await chat.sendMessage({
          message: message.content
        });
      } else {
        // For user messages
        await chat.sendMessage({
          message: message.content
        });
      }
    }

    // Create a stream
    const encoder = new TextEncoder();
    const stream = new ReadableStream({
      async start(controller) {
        try {
          // Generate a streaming response
          const response = await chat.sendMessageStream({
            message: messages[messages.length - 1].content,
            config: {
              temperature,
              maxOutputTokens: maxOutputTokens || 4096, // Use provided max tokens or default
            }
          });

          // Process each chunk
          for await (const chunk of response) {
            const text = chunk.text;
            if (text) {
              // Encode the chunk as a data event in the SSE format
              const data = `data: ${JSON.stringify({ choices: [{ delta: { content: text } }] })}\n\n`;
              controller.enqueue(encoder.encode(data));
            }
          }

          // Signal the end of the stream
          controller.enqueue(encoder.encode('data: [DONE]\n\n'));
          controller.close();
        } catch (error) {
          console.error('Error in Gemini streaming:', error);
          controller.error(error);
        }
      }
    });

    return stream;
  } catch (error) {
    console.error('Error generating Gemini streaming response:', error);
    throw error;
  }
}
