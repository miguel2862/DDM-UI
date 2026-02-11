import { phenomena } from '../../data/phenomena';
import { PhenomenonCard } from './PhenomenonCard';
import { useI18n } from '../../i18n';

interface PhenomenonGalleryProps {
  onLoad: (id: string) => void;
  loadingId?: string | null;
}

export function PhenomenonGallery({ onLoad, loadingId }: PhenomenonGalleryProps) {
  const { t } = useI18n();

  return (
    <div>
      <div className="flex items-center gap-3 mb-6">
        <h2 className="text-xl font-bold text-navy-100">{t.gallery.title}</h2>
        <span className="text-xs text-navy-500 font-medium px-2 py-0.5 rounded-full bg-navy-800/60 border border-navy-700/50">
          {phenomena.filter(p => p.available).length} / {phenomena.length} {t.gallery.available}
        </span>
      </div>
      <div className="grid grid-cols-1 md:grid-cols-2 xl:grid-cols-3 gap-5">
        {phenomena.map((p, i) => (
          <PhenomenonCard
            key={p.id}
            phenomenon={p}
            onLoad={onLoad}
            delay={i * 0.08}
            isLoading={loadingId === p.id}
          />
        ))}
      </div>
    </div>
  );
}
