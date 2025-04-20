import { getProgramById } from './config';
import { Program } from './types';
import { getEmbeddedContext as getEmbeddedContextFromModule } from './embedded-context';

// Cache for storing loaded context to avoid repeated lookups
const combinedContextCache: Record<string, string> = {};

/**
 * Load context for a program using the embedded context content
 * @param contextFiles Array of context file paths or a single context file path
 * @returns The content of the combined context
 */
export async function loadContext(contextFiles: string[] | string): Promise<string> {
  // Convert single file to array for consistent handling
  const filesArray = typeof contextFiles === 'string' ? [contextFiles] : contextFiles;

  // Generate a cache key based on the context files
  const cacheKey = filesArray.join('|');

  // Check if the combined context is already cached
  if (combinedContextCache[cacheKey]) {
    return combinedContextCache[cacheKey];
  }

  // Find the program that matches these context files
  const program = Object.values(getProgramById() as Record<string, Program>)
    .find(p => p.contextFiles &&
      JSON.stringify(p.contextFiles.sort()) === JSON.stringify(filesArray.sort()));

  if (!program) {
    console.error(`No program found for context files: ${filesArray.join(', ')}`);
    return `Context could not be loaded. No matching program found for the specified context files.`;
  }

  console.log(`Loading embedded context for program: ${program.id}`);

  // Get the embedded context for this program
  const contextContent = getEmbeddedContextFromModule(program.id);

  if (!contextContent) {
    console.error(`No embedded context found for program: ${program.id}`);
    return `Context could not be loaded. No embedded context found for program: ${program.id}.`;
  }

  console.log(`Successfully loaded embedded context for ${program.id} (${contextContent.length} bytes)`);

  // Cache the content
  combinedContextCache[cacheKey] = contextContent;
  return contextContent;
}

/**
 * Get the embedded context for a program directly
 * This is a convenience function that wraps the imported getEmbeddedContext
 * @param programId The ID of the program
 * @returns The embedded context content or undefined if not found
 */
export function getEmbeddedContext(programId: string): string | undefined {
  // Re-export from embedded-context.ts
  return getEmbeddedContextFromModule(programId);
}

/**
 * Get the system prompt with context for a specific program
 * @param programId The ID of the program
 * @param context The context content
 * @returns The system prompt with context
 */
export function getSystemPromptWithContext(programId: string, context: string): string {
  // Get the program configuration and global config
  const program = getProgramById(programId) as Program;
  const { loadConfig } = require('./config');
  const config = loadConfig();

  // Get program name and uppercase ID
  const programName = program?.name || programId.toUpperCase();
  const programIdUpper = programId.toUpperCase();

  // Get the system prompt template and replace placeholders
  // Load the system prompt from config.json if it's not available in the config object
  const rawConfig = require('../../config.json');
  let systemPromptTemplate = (config.systemPrompt || rawConfig.systemPrompt || '')
    .replace(/\{PROGRAM_ID\}/g, programIdUpper)
    .replace(/\{PROGRAM_NAME\}/g, programName);

  // Add program-specific extra system prompt if available
  if (program?.extraSystemPrompt) {
    systemPromptTemplate += `\n\n${program.extraSystemPrompt}`;
  }

  // Add additional context if available
  if (config.additionalContext) {
    systemPromptTemplate += `\n\nUSER CONTEXT:\n${config.additionalContext}`;
  }

  // Append the documentation to the system prompt
  return `${systemPromptTemplate}\n\nDOCUMENTATION:\n${context}`;
}
