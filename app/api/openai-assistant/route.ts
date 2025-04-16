import { NextRequest } from 'next/server';
import OpenAI from 'openai';

// this enables Edge Functions in Vercel
export const runtime = "edge";

// post a new message and stream OpenAI Assistant response
export async function POST(request: NextRequest) {
    try {
        // parse message from post
        const newMessage = await request.json();

        // create OpenAI client
        const openai = new OpenAI();

        // if no thread id then create a new openai thread
        let threadId = newMessage.threadId;
        if (!threadId) {
            const thread = await openai.beta.threads.create();
            threadId = thread.id;
        }

        // add new message to thread
        await openai.beta.threads.messages.create(
            threadId,
            {
                role: "user",
                content: newMessage.content
            }
        );

        // create a run
        const stream = await openai.beta.threads.runs.create(
            threadId, 
            {assistant_id: newMessage.assistantId, stream: true}
        );
        
        return new Response(stream.toReadableStream());
    } catch (error) {
        console.error('Error in OpenAI Assistant API:', error);
        return new Response(
            JSON.stringify({ error: 'An error occurred processing your request' }),
            { status: 500 }
        );
    }
}
