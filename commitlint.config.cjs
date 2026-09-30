module.exports = {
  extends: ["@commitlint/config-conventional"],
  rules: {
    "scope-enum": [
      2,
      "always",
      ["server", "desktop", "electron", "tooling", "ci", "docker", "deps"],
    ],
  },
}
