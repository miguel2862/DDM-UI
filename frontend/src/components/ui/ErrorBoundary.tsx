import { Component, type ReactNode } from 'react';
import { AlertTriangle, RotateCcw, Trash2 } from 'lucide-react';
import { useI18n } from '../../i18n';

interface Props {
  children: ReactNode;
}

interface State {
  hasError: boolean;
  error: Error | null;
  showDetails: boolean;
}

export class ErrorBoundary extends Component<Props, State> {
  constructor(props: Props) {
    super(props);
    this.state = { hasError: false, error: null, showDetails: false };
  }

  static getDerivedStateFromError(error: Error): Partial<State> {
    return { hasError: true, error };
  }

  componentDidCatch(error: Error, info: { componentStack?: string | null }) {
    console.error('[DDM-UI] Uncaught error:', error, info.componentStack);
  }

  handleReload = () => {
    window.location.reload();
  };

  handleReset = () => {
    try {
      localStorage.removeItem('ddm-ui-autosave');
    } catch { /* ignore */ }
    window.location.reload();
  };

  render() {
    if (this.state.hasError) {
      const t = useI18n.getState().t;
      return (
        <div className="min-h-screen bg-gradient-to-b from-slate-900 via-slate-800 to-slate-900 flex items-center justify-center p-8">
          <div className="max-w-md w-full space-y-6 text-center">
            <div className="mx-auto w-16 h-16 rounded-2xl bg-rose-500/20 flex items-center justify-center">
              <AlertTriangle size={32} className="text-rose-400" />
            </div>

            <div>
              <h1 className="text-2xl font-bold text-white mb-2">
                {t.errorBoundary.title}
              </h1>
              <p className="text-sm text-slate-400">
                {t.errorBoundary.description}
              </p>
            </div>

            <div className="flex flex-col gap-3">
              <button
                onClick={this.handleReload}
                className="w-full flex items-center justify-center gap-2 px-6 py-3 rounded-xl bg-cyan-500 text-white font-semibold text-sm hover:bg-cyan-400 transition-colors"
              >
                <RotateCcw size={16} />
                {t.errorBoundary.reload}
              </button>
              <button
                onClick={this.handleReset}
                className="w-full flex items-center justify-center gap-2 px-6 py-3 rounded-xl bg-slate-700 text-slate-300 font-semibold text-sm hover:bg-slate-600 transition-colors"
              >
                <Trash2 size={16} />
                {t.errorBoundary.reset}
              </button>
            </div>

            {this.state.error && (
              <div>
                <button
                  onClick={() => this.setState(s => ({ showDetails: !s.showDetails }))}
                  className="text-xs text-slate-500 hover:text-slate-400 transition-colors"
                >
                  {this.state.showDetails ? '▾' : '▸'} {t.errorBoundary.details}
                </button>
                {this.state.showDetails && (
                  <pre className="mt-2 p-3 rounded-lg bg-slate-800 text-left text-[10px] text-rose-300 overflow-auto max-h-40 border border-slate-700">
                    {this.state.error.message}
                    {'\n\n'}
                    {this.state.error.stack}
                  </pre>
                )}
              </div>
            )}
          </div>
        </div>
      );
    }

    return this.props.children;
  }
}
