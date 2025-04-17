'use client'

import { useState, useRef, useEffect } from "react";
import { AiOutlineUser, AiOutlineRobot, AiOutlineSend } from "react-icons/ai";
import Markdown from 'react-markdown';
import rehypeHighlight from 'rehype-highlight';
import rehypeRaw from 'rehype-raw';
import rehypeKatex from 'rehype-katex';
import remarkGfm from 'remark-gfm';
import remarkMath from 'remark-math';
import { Message } from '@/app/utils/types';
// We'll handle syntax highlighting styles in globals.css
import CopyButton from './copy-button';

interface ChatProps {
  programId: string;
  greeting: string;
}

export default function Chat({
  programId,
  greeting = "How can I help you?",
}: ChatProps) {
  const [isLoading, setIsLoading] = useState(false);
  const [prompt, setPrompt] = useState("");
  const [messages, setMessages] = useState<Message[]>([]);
  const [streamingContent, setStreamingContent] = useState("");
  const messageId = useRef(0);

  // Create a ref for the greeting message to avoid hydration mismatches
  const greetingMessageRef = useRef<Message>({
    id: "greeting",
    role: "assistant",
    content: `Welcome to the ${programId} help assistant. ${greeting}`,
    createdAt: new Date(0), // Use a consistent date initially
  });

  // Update the greeting message when the program changes
  useEffect(() => {
    console.log('Program changed to:', programId);
    greetingMessageRef.current = {
      id: "greeting",
      role: "assistant",
      content: `Welcome to the ${programId} help assistant. ${greeting}`,
      createdAt: new Date()
    };
    // Reset messages when program changes
    setMessages([]);
  }, [programId, greeting]);

  async function handleSubmit(e: React.FormEvent<HTMLFormElement>) {
    e.preventDefault();

    // Clear streaming content
    setStreamingContent("");

    // Add busy indicator
    setIsLoading(true);

    // Add user message to list of messages
    messageId.current++;
    // Create a timestamp once for this message
    const messageTimestamp = new Date();

    const userMessage = {
      id: messageId.current.toString(),
      role: "user",
      content: prompt,
      createdAt: messageTimestamp,
    };

    setMessages([...messages, userMessage]);
    setPrompt("");

    // Prepare messages for API
    const apiMessages = [
      ...messages.map(msg => ({
        role: msg.role,
        content: msg.content
      })),
      {
        role: userMessage.role,
        content: userMessage.content
      }
    ];

    try {
      console.log('Sending request to API with programId:', programId);

      // Use the non-streaming API for testing
      const response = await fetch('/api/chat-simple', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          programId,
          messages: apiMessages,
        }),
      });

      console.log('API response status:', response.status);

      // Check if the response is ok
      if (!response.ok) {
        const errorData = await response.json().catch(() => ({}));
        const errorDetails = errorData.details || '';
        const errorMessage = errorData.error || response.statusText;

        // Provide more specific error messages based on the error type
        if (errorMessage.includes('API key')) {
          throw new Error(`OpenAI API key error: ${errorDetails}`);
        } else if (errorMessage.includes('model')) {
          throw new Error(`Model configuration error: ${errorDetails}`);
        } else {
          throw new Error(`API responded with status ${response.status}: ${errorMessage}`);
        }
      }

      // For non-streaming API, parse the JSON response
      const responseData = await response.json();
      console.log('Response data:', responseData);

      // Extract the content from the response
      const finalContent = responseData.content || "No content received from API";
      console.log('Final content:', finalContent);

      // Add assistant message to list of messages
      messageId.current++;

      // Use the same timestamp for consistency
      setMessages(prevMessages => [
        ...prevMessages,
        {
          id: messageId.current.toString(),
          role: "assistant",
          content: finalContent,
          createdAt: messageTimestamp, // Reuse the same timestamp
        }
      ]);
    } catch (error) {
      console.error('Error in chat:', error);

      // Display error message to the user
      messageId.current++;
      const errorMessage = error instanceof Error ? error.message : 'Unknown error';
      // Create a timestamp for the error message
      const errorTimestamp = new Date();
      setMessages(prevMessages => [
        ...prevMessages,
        {
          id: messageId.current.toString(),
          role: "assistant",
          content: `I'm sorry, there was an error processing your request: ${errorMessage}. Please try again later.`,
          createdAt: errorTimestamp,
        }
      ]);
    } finally {
      // Remove busy indicator
      setIsLoading(false);
      setStreamingContent("");
    }
  }

  // Handles changes to the prompt input field
  function handlePromptChange(e: React.ChangeEvent<HTMLInputElement>) {
    setPrompt(e.target.value);
  }

  return (
    <div className="flex flex-col bg-white dark:bg-gray-800">
      <div className="p-2">
        <ChatMessage
          message={greetingMessageRef.current}
        />
        {messages.map(m =>
          <ChatMessage
            key={m.id}
            message={m}
          />
        )}
        {isLoading && streamingContent && (
          <ChatMessage
            message={{
              id: "streaming",
              role: "assistant",
              content: streamingContent,
              createdAt: new Date(0), // Use a consistent date to avoid hydration issues
            }}
          />
        )}
      </div>
      <form onSubmit={handleSubmit} className="p-2 border-t border-gray-300 dark:border-gray-700 flex">
        <input
          disabled={isLoading}
          autoFocus
          className="border border-gray-300 dark:border-gray-700 rounded w-full py-2 px-3 text-gray-700 dark:text-gray-300 bg-white dark:bg-gray-800"
          onChange={handlePromptChange}
          value={prompt}
          placeholder={isLoading ? "Thinking..." : "Ask a question..."} />
        {isLoading ?
          <button
            disabled
            className="ml-2 bg-blue-500 text-white font-bold py-2 px-4 rounded focus:outline-none">
            <ChatSpinner />
          </button>
          :
          <button
            disabled={prompt.length === 0}
            className="ml-2 bg-blue-500 hover:bg-blue-700 text-white font-bold py-2 px-4 rounded focus:outline-none">
            <AiOutlineSend />
          </button>
        }
      </form>
    </div>
  );
}

function ChatMessage({ message }: { message: Message }) {
  function displayRole(roleName: string) {
    switch (roleName) {
      case "user":
        return <AiOutlineUser className="text-gray-600 dark:text-gray-300" />;
      case "assistant":
        return <AiOutlineRobot className="text-gray-600 dark:text-gray-300" />;
      default:
        return null;
    }
  }

  const isUser = message.role === "user";

  return (
    <div className={`flex rounded-lg text-gray-700 dark:text-gray-200 px-4 py-3 my-2 shadow-sm border ${isUser ? 'bg-blue-50 border-blue-100 dark:bg-blue-900 dark:border-blue-800' : 'bg-white border-gray-200 dark:bg-gray-800 dark:border-gray-700'}`}>
      <div className="text-3xl flex-shrink-0 flex items-start pt-1">
        {displayRole(message.role)}
      </div>
      <div className="mx-3 text-left overflow-auto w-full prose dark:prose-invert max-w-none">
        <Markdown
          remarkPlugins={[remarkGfm, remarkMath]}
          rehypePlugins={[rehypeRaw, rehypeHighlight, rehypeKatex]}
          components={{
            code(props) {
              const { node, inline, className, children } = props as {
                node?: any;
                inline?: boolean;
                className?: string;
                children: React.ReactNode;
              };
              // We don't need to use match, so we can remove it
              // const match = /language-(\w+)/.exec(className || '');

              // Extract text for copy button
              const extractText = (node: any): string => {
                if (typeof node === 'string') return node;
                if (Array.isArray(node)) return node.map(extractText).join('');
                if (node && typeof node === 'object' && 'props' in node) {
                  return extractText(node.props?.children || '');
                }
                return '';
              };

              // For inline code, keep it simple
              if (inline) {
                return <code className={`${className} inline-code`} {...props}>{children}</code>;
              }

              // For code blocks, check if it's inside a pre tag
              // This is the most reliable way to identify actual code blocks vs standalone code tags
              const isInPreTag = node &&
                'parentNode' in node &&
                node.parentNode &&
                'tagName' in node.parentNode &&
                node.parentNode.tagName &&
                node.parentNode.tagName.toLowerCase() === 'pre';

              if (isInPreTag) {
                // For actual code blocks (inside pre tags), add copy button
                return (
                  <>
                    <code className={className} {...props}>{children}</code>
                    <CopyButton text={extractText(children)} />
                  </>
                );
              } else {
                // For standalone code tags that aren't inline (rare case), just render the code
                return <code className={className} {...props}>{children}</code>;
              }
            }
          }}
        >
          {message.content}
        </Markdown>
      </div>
    </div>
  );
}

function ChatSpinner() {
  return (
    <svg className="animate-spin h-5 w-5" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24">
      <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4"></circle>
      <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path>
    </svg>
  );
}
