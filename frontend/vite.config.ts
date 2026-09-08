import { defineConfig, loadEnv } from 'vite'
import react from '@vitejs/plugin-react'
import tailwindcss from '@tailwindcss/vite'

export default defineConfig(({ mode }) => {
  const env = loadEnv(mode, process.cwd(), '')
  return {
    plugins: [react(), tailwindcss()],
    // Relative paths for Electron file:// protocol compatibility
    base: './',
    server: {
      port: 5173,
      proxy: {
        '/api': {
          target: env.VITE_API_TARGET || 'http://localhost:8000',
          changeOrigin: true,
        },
      },
    },
    build: {
      rollupOptions: {
        output: {
          manualChunks(id) {
            if (!id.includes('node_modules')) return undefined
            if (id.includes('@xyflow') || id.includes('d3-')) return 'network-vendor'
            if (id.includes('recharts') || id.includes('victory-vendor')) return 'charts-vendor'
            if (id.includes('framer-motion') || id.includes('motion-')) return 'motion-vendor'
            if (id.includes('react') || id.includes('zustand') || id.includes('redux')) return 'react-vendor'
            return undefined
          },
        },
      },
    },
  }
})
