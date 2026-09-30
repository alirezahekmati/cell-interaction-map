import { defineConfig } from 'vite';
// relative base so the site works under https://<user>.github.io/<repo>/
export default defineConfig({ base: './', build: { outDir: '../docs', emptyOutDir: true } });
