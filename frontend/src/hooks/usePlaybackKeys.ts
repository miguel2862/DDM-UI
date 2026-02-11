import { useEffect } from 'react';
import { useSimStore } from '../stores/useSimStore';

/**
 * Global keyboard shortcuts for simulation playback:
 * - Space: Play/Pause
 * - Left Arrow: Step backward
 * - Right Arrow: Step forward
 * - Shift+Right: Skip +10 steps
 * - Home: Go to start
 * - End: Go to end
 * - 1-4: Set speed (1x, 2x, 5x, 10x)
 *
 * Only active on /simulation and /results pages.
 */
export function usePlaybackKeys() {
  useEffect(() => {
    const handler = (e: KeyboardEvent) => {
      const state = useSimStore.getState();

      // Only activate when simulation is complete and we have data
      if (state.simStatus !== 'complete' || !state.simulationResults?.length) return;

      // Only activate on simulation and results pages (HashRouter uses hash)
      const hash = window.location.hash;
      if (!hash.includes('/simulation') && !hash.includes('/results')) return;

      // Don't capture keys when typing in inputs or contentEditable
      const el = e.target as HTMLElement;
      const tag = el.tagName;
      if (tag === 'INPUT' || tag === 'TEXTAREA' || tag === 'SELECT' || el.isContentEditable) return;

      const maxIndex = (state.simulationResults[0]?.length || 1) - 1;

      switch (e.key) {
        case ' ':
          e.preventDefault();
          state.setIsPlaying(!state.isPlaying);
          break;
        case 'ArrowLeft':
          e.preventDefault();
          state.setIsPlaying(false);
          state.setPlaybackIndex(Math.max(0, state.playbackIndex - 1));
          break;
        case 'ArrowRight':
          e.preventDefault();
          state.setIsPlaying(false);
          if (e.shiftKey) {
            state.setPlaybackIndex(Math.min(maxIndex, state.playbackIndex + 10));
          } else {
            state.setPlaybackIndex(Math.min(maxIndex, state.playbackIndex + 1));
          }
          break;
        case 'Home':
          e.preventDefault();
          state.setIsPlaying(false);
          state.setPlaybackIndex(0);
          break;
        case 'End':
          e.preventDefault();
          state.setIsPlaying(false);
          state.setPlaybackIndex(maxIndex);
          break;
        case '1':
          state.setPlaybackSpeed(1);
          break;
        case '2':
          state.setPlaybackSpeed(2);
          break;
        case '3':
          state.setPlaybackSpeed(5);
          break;
        case '4':
          state.setPlaybackSpeed(10);
          break;
      }
    };

    window.addEventListener('keydown', handler);
    return () => window.removeEventListener('keydown', handler);
  }, []);
}
