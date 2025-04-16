import { NextRequest } from 'next/server';
import { getEmbeddedContext } from '@/app/utils/context';

export async function GET(_request: NextRequest, context: { params: Promise<{ programId: string }> }) {
  try {
    // Get the programId from params
    const { programId } = await context.params;

    // Get the embedded context for the specified program
    const contextContent = getEmbeddedContext(programId);

    if (!contextContent) {
      return new Response(
        JSON.stringify({
          error: `No embedded context found for program: ${programId}`,
          timestamp: new Date().toISOString()
        }),
        {
          status: 404,
          headers: {
            'Content-Type': 'application/json'
          }
        }
      );
    }

    // Return the context as plain text
    return new Response(contextContent, {
      status: 200,
      headers: {
        'Content-Type': 'text/markdown'
      }
    });
  } catch (error) {
    console.error(`Error accessing embedded context:`, error);

    return new Response(
      JSON.stringify({
        error: 'An error occurred accessing embedded context',
        details: error instanceof Error ? error.message : String(error),
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
