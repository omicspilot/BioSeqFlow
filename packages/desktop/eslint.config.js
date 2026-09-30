import js from "@eslint/js"
import tseslint from "typescript-eslint"
import reactHooks from "eslint-plugin-react-hooks"
import jsxA11y from "eslint-plugin-jsx-a11y"
import prettierConfig from "eslint-config-prettier"

export default tseslint.config(
  // 1. Global ignores — never lint generated or compiled output
  { ignores: ["dist/", "out/", "node_modules/", "src/api/"] },

  // 2. Base JS rules for all files
  js.configs.recommended,

  // 3. TypeScript strict rules, scoped to .ts/.tsx files only
  {
    files: ["**/*.ts", "**/*.tsx"],
    extends: tseslint.configs.strict,
    languageOptions: {
      parserOptions: {
        project: true,
        tsconfigRootDir: import.meta.dirname,
      },
    },
    plugins: {
      "react-hooks": reactHooks,
    },
    rules: {
      ...reactHooks.configs.recommended.rules,
      "no-console": "error",
    },
  },

  // 4. Accessibility rules, scoped to .tsx only (JSX lives there)
  {
    files: ["**/*.tsx"],
    ...jsxA11y.flatConfigs.recommended,
  },

  // 5. Prettier must be last — disables any ESLint rules that conflict with formatting
  prettierConfig,
)
