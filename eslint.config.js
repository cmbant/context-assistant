const nextPlugin = require('@next/eslint-plugin-next');

module.exports = [
  {
    ignores: ['.next/**', 'node_modules/**', 'dist/**', 'build/**', '*.config.js'],
  },
  {
    plugins: {
      '@next/next': nextPlugin,
    },
    rules: {
      ...nextPlugin.configs['core-web-vitals'].rules,
    },
  },
];
