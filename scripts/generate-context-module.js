const fs = require('fs');
const path = require('path');

// Read programs from config.json
const configPath = path.join(__dirname, '..', 'config.json');
const config = JSON.parse(fs.readFileSync(configPath, 'utf8'));
const programs = config.programs;

console.log('Generating embedded context module...');

// Source directory for context files
const contextDir = path.join(__dirname, '..', 'context');

// Generate combined context content
const combinedContexts = {};

for (const program of programs) {
  try {
    const contents = program.contextFiles.map(file => {
      const filePath = path.join(contextDir, file);
      if (!fs.existsSync(filePath)) {
        console.error(`Context file not found: ${filePath}`);
        return `Context file ${file} not found.`;
      }
      return fs.readFileSync(filePath, 'utf8');
    });

    const combinedContent = contents.join('\n\n---\n\n');
    combinedContexts[program.id] = combinedContent;

    console.log(`Generated combined context for ${program.id} (${combinedContent.length} bytes)`);
  } catch (error) {
    console.error(`Error generating combined context for ${program.id}:`, error);
  }
}

// Generate TypeScript module with embedded context
const tsContent = `// This file is auto-generated. Do not edit directly.
// Generated on ${new Date().toISOString()}

/**
 * Embedded context content for each program
 * This avoids having to load context files at runtime
 */
export const embeddedContexts: Record<string, string> = {
${Object.entries(combinedContexts).map(([id, content]) => {
  // Use JSON.stringify to properly escape all special characters
  return `  '${id}': ${JSON.stringify(content)},`;
}).join('\n')}
};

/**
 * Get the embedded context for a program
 * @param programId The ID of the program
 * @returns The embedded context content or undefined if not found
 */
export function getEmbeddedContext(programId: string): string | undefined {
  return embeddedContexts[programId];
}
`;

// Write the TypeScript file
const outputPath = path.join(__dirname, '..', 'app', 'utils', 'embedded-context.ts');
fs.writeFileSync(outputPath, tsContent, 'utf8');

console.log(`Generated embedded context module at ${outputPath}`);
