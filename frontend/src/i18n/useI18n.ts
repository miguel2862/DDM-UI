import { create } from 'zustand';
import { persist } from 'zustand/middleware';
import { en } from './en';
import { es } from './es';
import type { Translations } from './en';

export type Language = 'en' | 'es';

interface I18nStore {
  language: Language;
  setLanguage: (lang: Language) => void;
  t: Translations;
}

const translations: Record<Language, Translations> = { en, es };

export const useI18n = create<I18nStore>()(
  persist(
    (set) => ({
      language: 'en',
      t: en,
      setLanguage: (lang: Language) =>
        set({ language: lang, t: translations[lang] }),
    }),
    {
      name: 'ddm-ui-language',
      partialize: (state) => ({ language: state.language }),
      onRehydrateStorage: () => {
        return (state?: I18nStore) => {
          if (state) {
            state.t = translations[state.language];
          }
        };
      },
    }
  )
);
