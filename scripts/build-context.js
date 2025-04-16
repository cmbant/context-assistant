const fs = require('fs');
const path = require('path');

// Read programs from config.json
const configPath = path.join(__dirname, '..', 'config.json');
const config = JSON.parse(fs.readFileSync(configPath, 'utf8'));
const programs = config.programs;

console.log('Using programs:', programs);

// Ensure public/context directory exists
const publicContextDir = path.join(__dirname, '..', 'public', 'context');
if (!fs.existsSync(publicContextDir)) {
  fs.mkdirSync(publicContextDir, { recursive: true });
}

// Copy individual context files to public/context
const contextDir = path.join(__dirname, '..', 'context');
const contextFiles = fs.readdirSync(contextDir);

for (const file of contextFiles) {
  const sourcePath = path.join(contextDir, file);
  const destPath = path.join(publicContextDir, file);

  if (fs.statSync(sourcePath).isFile()) {
    fs.copyFileSync(sourcePath, destPath);
    console.log(`Copied ${file} to public/context`);
  }
}

// Generate combined context files
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
    const outputPath = path.join(publicContextDir, program.combinedContextFile);

    fs.writeFileSync(outputPath, combinedContent, 'utf8');
    console.log(`Generated combined context file for ${program.id}: ${program.combinedContextFile}`);
  } catch (error) {
    console.error(`Error generating combined context for ${program.id}:`, error);
  }
}

console.log('Context files built successfully');
