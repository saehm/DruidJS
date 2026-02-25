import { playwright } from "@vitest/browser-playwright";
import { defineConfig } from "vitest/config";

export default defineConfig({
    test: {
        projects: [
            {
                test: {
                    name: "node",
                    environment: "node",
                    typecheck: { enabled: true },
                    include: ["test/**/*.js", "test/**/*.test.js"],
                    exclude: ["test/test_index.js", "test/utils/**"],
                    coverage: { provider: "v8" },
                },
            },
            {
                test: {
                    name: "browser-tests",
                    include: ["test/**/*.js", "test/**/*.test.js"],
                    exclude: ["test/test_index.js", "test/utils/**"],
                    browser: {
                        enabled: true,
                        provider: playwright(),
                        instances: [{ browser: "chromium" }, { browser: "firefox" }, { browser: "webkit" }],
                        headless: true,
                    },
                },
            },
        ],
    },
});
