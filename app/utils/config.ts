import { Config, Program } from './types';

// Import the config directly
import rawConfigJson from '../../config.json';

// Add type assertion to include the new greeting property
const configJson = rawConfigJson as unknown as Partial<Config>;

// Default configuration for non-program settings
const defaultConfig: Partial<Config> = {
  showContextLink: true, // Set to false to hide the context link
  simpleMode: false, // Set to true to hide top panels when only one code and model
  greeting: "How can I help you?", // Default greeting message
  apiKeys: {
    // These are just placeholders, actual keys should be set in environment variables
    openai: process.env.OPENAI_API_KEY,
    gemini: process.env.GEMINI_API_KEY
  }
  // availableModels and defaultModelId are now loaded from config.json
};

/**
 * Load the configuration
 * @returns The configuration object
 */
export function loadConfig(): Config {
  // Merge with default config to ensure all fields are present
  return {
    ...defaultConfig,
    programs: configJson.programs,
    defaultProgram: configJson.defaultProgram,
    showContextLink: configJson.showContextLink !== undefined ? configJson.showContextLink : defaultConfig.showContextLink,
    simpleMode: configJson.simpleMode !== undefined ? configJson.simpleMode : defaultConfig.simpleMode,
    greeting: configJson.greeting !== undefined ? configJson.greeting : defaultConfig.greeting,
    defaultModelId: configJson.defaultModelId || 'gemini/gemini-2.0-flash',
    availableModels: configJson.availableModels || [],
    // Load the new flags, defaulting to true if not present in config.json
    useDirectOpenAIKey: configJson.useDirectOpenAIKey ?? true,
    useDirectGeminiKey: configJson.useDirectGeminiKey ?? true,
    apiKeys: defaultConfig.apiKeys || {}
  } as Config;
}

/**
 * Get a program by its ID or all programs if no ID is provided
 * @param programId The ID of the program (optional)
 * @returns The program object, all programs as a record, or undefined if not found
 */
export function getProgramById(programId?: string): Program | Record<string, Program> | undefined {
  const config = loadConfig();

  // If no programId is provided, return all programs as a record
  if (!programId) {
    const programsRecord: Record<string, Program> = {};
    config.programs.forEach(program => {
      programsRecord[program.id] = program;
    });
    return programsRecord;
  }

  // Otherwise, return the specific program
  return config.programs.find(program => program.id === programId);
}

/**
 * Get the default program
 * @returns The default program object
 */
export function getDefaultProgram(): Program {
  const config = loadConfig();
  const defaultProgram = config.programs.find(program => program.id === config.defaultProgram);
  return defaultProgram || config.programs[0];
}

/**
 * Parse a model ID in the format 'provider/model-name'
 * @param modelId The model ID to parse
 * @returns An object with provider and modelName properties
 */
export function parseModelId(modelId: string): { provider: string; modelName: string } {
  const parts = modelId.split('/');
  if (parts.length !== 2) {
    throw new Error(`Invalid model ID format: ${modelId}. Expected format: 'provider/model-name'`);
  }
  return {
    provider: parts[0],
    modelName: parts[1]
  };
}


