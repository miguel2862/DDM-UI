import { HelpCircle, BookOpen, ExternalLink } from 'lucide-react';
import { Card } from '../components/ui/Card';
import { PageTransition } from '../components/ui/Skeleton';
import { useI18n } from '../i18n';

export function Help() {
  const { t } = useI18n();

  return (
    <PageTransition>
    <div className="max-w-4xl mx-auto space-y-8">
      {/* Page header */}
      <div className="flex items-center gap-3">
        <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-cyan-500 to-teal-500 flex items-center justify-center">
          <HelpCircle size={20} className="text-white" />
        </div>
        <div>
          <h1 className="text-2xl font-bold text-slate-800">{t.help.pageTitle}</h1>
          <p className="text-sm text-slate-500">{t.help.pageSubtitle}</p>
        </div>
      </div>

      {/* About the Model */}
      <Card className="p-6">
        <div className="flex items-center gap-2 mb-4">
          <BookOpen size={18} className="text-cyan-600" />
          <h2 className="text-lg font-bold text-slate-800">{t.help.aboutTitle}</h2>
        </div>
        <div className="space-y-3 text-sm text-slate-600 leading-relaxed">
          <p>{t.help.aboutP1}</p>
          <p>{t.help.aboutP2}</p>
          <p>{t.help.aboutP3}</p>
        </div>
      </Card>




      {/* Network Architecture */}
      <Card className="p-6">
        <div className="flex items-center gap-2 mb-4">
          <BookOpen size={18} className="text-cyan-600" />
          <h2 className="text-lg font-bold text-slate-800">{t.help.networkTitle}</h2>
        </div>
        <div className="space-y-3 text-sm text-slate-600 leading-relaxed">
          <p>{t.help.networkIntro}</p>
          <ul className="space-y-2 ml-4">
            <li>
              <span className="font-semibold text-slate-700">{t.help.layerUS}</span> &mdash;
              {t.help.layerUSDesc}
            </li>
            <li>
              <span className="font-semibold text-slate-700">{t.help.layerPS}</span> &mdash;
              {t.help.layerPSDesc}
            </li>
            <li>
              <span className="font-semibold text-slate-700">{t.help.layerAS}</span> &mdash;
              {t.help.layerASDesc}
            </li>
            <li>
              <span className="font-semibold text-slate-700">{t.help.layerH}</span> &mdash;
              {t.help.layerHDesc}
            </li>
            <li>
              <span className="font-semibold text-slate-700">{t.help.layerAM}</span> &mdash;
              {t.help.layerAMDesc}
            </li>
            <li>
              <span className="font-semibold text-slate-700">{t.help.layerPM}</span> &mdash;
              {t.help.layerPMDesc}
            </li>
            <li>
              <span className="font-semibold text-slate-700">{t.help.layerD}</span> &mdash;
              {t.help.layerDDesc}
            </li>
          </ul>
        </div>
      </Card>

      {/* Naming Convention */}
      <Card className="p-6">
        <div className="flex items-center gap-2 mb-4">
          <BookOpen size={18} className="text-cyan-600" />
          <h2 className="text-lg font-bold text-slate-800">{t.help.namingTitle}</h2>
        </div>
        <div className="space-y-3 text-sm text-slate-600 leading-relaxed">
          <p>{t.help.namingIntro}</p>
          <ul className="space-y-2 ml-4">
            <li>
              <span className="font-mono font-semibold text-slate-700">S..1</span> &mdash;{' '}
              {t.help.namingSdd}
            </li>
            <li>
              <span className="font-mono font-semibold text-slate-700">M..1</span> &mdash;{' '}
              {t.help.namingMdd}
            </li>
            <li>
              <span className="font-mono font-semibold text-slate-700">M.1</span> &mdash;{' '}
              {t.help.namingMd}
            </li>
            <li>
              <span className="font-mono font-semibold text-slate-700">H1, H2</span> &mdash;{' '}
              {t.help.namingH}
            </li>
            <li>
              <span className="font-mono font-semibold text-slate-700">D</span> &mdash;{' '}
              {t.help.namingD}
            </li>
            <li>
              <span className="font-mono font-semibold text-slate-700">US</span> &mdash;{' '}
              {t.help.namingUS}
            </li>
            <li>
              <span className="font-mono font-semibold text-slate-700">S1, S2, ...</span> &mdash;{' '}
              {t.help.namingS}
            </li>
          </ul>
        </div>
      </Card>

      {/* Parameters */}
      <Card className="p-6">
        <div className="flex items-center gap-2 mb-4">
          <BookOpen size={18} className="text-cyan-600" />
          <h2 className="text-lg font-bold text-slate-800">{t.help.freeParamsTitle}</h2>
        </div>
        <div className="space-y-3 text-sm text-slate-600 leading-relaxed">
          <p>{t.help.freeParamsIntro}</p>

          <h3 className="font-semibold text-slate-700 mt-4">{t.help.npeParams}</h3>
          <ul className="space-y-2 ml-4">
            <li>
              <span className="font-semibold text-slate-700">&mu; (mu)</span> &mdash;{' '}
              {t.help.paramMu}
            </li>
            <li>
              <span className="font-semibold text-slate-700">&sigma; (sigma)</span> &mdash;{' '}
              {t.help.paramSigma}
            </li>
            <li>
              <span className="font-semibold text-slate-700">&tau; (Temporal Summation)</span> &mdash;{' '}
              {t.help.paramTau}
            </li>
            <li>
              <span className="font-semibold text-slate-700">&kappa; (Activation Decay)</span> &mdash;{' '}
              {t.help.paramKappa}
            </li>
            <li>
              <span className="font-semibold text-slate-700">Logistic Slope</span> &mdash;{' '}
              {t.help.paramLogistic}
            </li>
          </ul>

          <h3 className="font-semibold text-slate-700 mt-4">{t.help.connParams}</h3>
          <ul className="space-y-2 ml-4">
            <li>
              <span className="font-semibold text-slate-700">Weight</span> &mdash;{' '}
              {t.help.paramWeight}
            </li>
            <li>
              <span className="font-semibold text-slate-700">&alpha; (alpha)</span> &mdash;{' '}
              {t.help.paramAlpha}
            </li>
            <li>
              <span className="font-semibold text-slate-700">&beta; (beta)</span> &mdash;{' '}
              {t.help.paramBeta}
            </li>
            <li>
              <span className="font-semibold text-slate-700">&alpha;&#x2032; (alpha prime)</span> &mdash;{' '}
              {t.help.paramAlphaPrime}
            </li>
            <li>
              <span className="font-semibold text-slate-700">&beta;&#x2032; (beta prime)</span> &mdash;{' '}
              {t.help.paramBetaPrime}
            </li>
          </ul>
        </div>
      </Card>

      {/* How to Use */}
      <Card className="p-6">
        <div className="flex items-center gap-2 mb-4">
          <BookOpen size={18} className="text-cyan-600" />
          <h2 className="text-lg font-bold text-slate-800">{t.help.howToUseTitle}</h2>
        </div>
        <div className="space-y-3 text-sm text-slate-600 leading-relaxed">
          <ol className="space-y-3 ml-4 list-decimal list-outside">
            <li>
              <span className="font-semibold text-slate-700">{t.help.step1Title}</span>{' '}
              {t.help.step1Desc}
            </li>
            <li>
              <span className="font-semibold text-slate-700">{t.help.step2Title}</span>{' '}
              {t.help.step2Desc}
            </li>
            <li>
              <span className="font-semibold text-slate-700">{t.help.step3Title}</span>{' '}
              {t.help.step3Desc}
            </li>
            <li>
              <span className="font-semibold text-slate-700">{t.help.step4Title}</span>{' '}
              {t.help.step4Desc}
            </li>
            <li>
              <span className="font-semibold text-slate-700">{t.help.step5Title}</span>{' '}
              {t.help.step5Desc}
            </li>
          </ol>
        </div>
      </Card>

      {/* References */}
      <Card className="p-6">
        <div className="flex items-center gap-2 mb-4">
          <BookOpen size={18} className="text-cyan-600" />
          <h2 className="text-lg font-bold text-slate-800">{t.help.referencesTitle}</h2>
        </div>
        <div className="space-y-4 text-sm text-slate-600 leading-relaxed">
          <div className="pl-6 -indent-6">
            <p>
              Aguayo-Mendoza, M., & Dos Santos, C.V. (2025). DDM-UI: A user interface in R for the
              discrepancy diffuse model in behavioral research.{' '}
              <em>Behavior Research Methods, 57</em>, 128.
            </p>
            <a href="https://doi.org/10.3758/s13428-025-02648-9" target="_blank" rel="noopener noreferrer"
              className="inline-flex items-center gap-1 text-cyan-600 hover:text-cyan-700 underline text-xs mt-1">
              <ExternalLink size={10} /> doi:10.3758/s13428-025-02648-9
            </a>
          </div>

          <div className="pl-6 -indent-6">
            <p>
              Donahoe, J.W., Burgos, J.E., & Palmer, D.C. (1993). A selectionist approach to
              reinforcement.{' '}
              <em>Journal of the Experimental Analysis of Behavior, 60</em>(1), 17&ndash;40.
            </p>
            <a href="https://doi.org/10.1901/jeab.1993.60-17" target="_blank" rel="noopener noreferrer"
              className="inline-flex items-center gap-1 text-cyan-600 hover:text-cyan-700 underline text-xs mt-1">
              <ExternalLink size={10} /> doi:10.1901/jeab.1993.60-17
            </a>
          </div>

          <div className="pl-6 -indent-6">
            <p>Aguayo-Mendoza, M., Buriticá, J., & Burgos, J.E. (2024). Autoshaped impulsivity: Some explorations with a neural network model. <em>Behavioural Processes, 218</em>, 105040.</p>
            <a href="https://doi.org/10.1016/j.beproc.2024.105040" target="_blank" rel="noopener noreferrer"
              className="inline-flex items-center gap-1 text-cyan-600 hover:text-cyan-700 underline text-xs mt-1">
              <ExternalLink size={10} /> doi:10.1016/j.beproc.2024.105040
            </a>
          </div>

          <div className="pl-6 -indent-6">
            <p>Lesaint, F., Sigaud, O., & Khamassi, M. (2014). Accounting for negative automaintenance in pigeons. <em>PLOS ONE, 9</em>, e111050.</p>
            <a href="https://doi.org/10.1371/journal.pone.0111050" target="_blank" rel="noopener noreferrer"
              className="inline-flex items-center gap-1 text-cyan-600 hover:text-cyan-700 underline text-xs mt-1">
              <ExternalLink size={10} /> doi:10.1371/journal.pone.0111050
            </a>
          </div>

          <div className="pl-6 -indent-6">
            <p>Jaskir, A., & Frank, M.J. (2023). On the normative advantages of dopamine and striatal opponency for learning and choice. <em>eLife, 12</em>, e85107.</p>
            <a href="https://doi.org/10.7554/eLife.85107" target="_blank" rel="noopener noreferrer"
              className="inline-flex items-center gap-1 text-cyan-600 hover:text-cyan-700 underline text-xs mt-1">
              <ExternalLink size={10} /> doi:10.7554/eLife.85107
            </a>
          </div>

        </div>
      </Card>

      {/* Studies using the model */}
      <Card className="p-6">
        <div className="flex items-center gap-2 mb-4">
          <BookOpen size={18} className="text-cyan-600" />
          <h2 className="text-lg font-bold text-slate-800">{t.help.studiesTitle}</h2>
        </div>
        <div className="space-y-2 text-sm text-slate-600 leading-relaxed">
          <p className="mb-3">
            {t.help.studiesIntro}
          </p>
          <div className="overflow-x-auto">
            <table className="w-full text-xs">
              <thead>
                <tr className="border-b border-slate-200">
                  <th className="text-left py-2 pr-3 font-semibold text-slate-700">{t.help.studyCol}</th>
                  <th className="text-left py-2 pr-3 font-semibold text-slate-700">{t.help.phenomenonCol}</th>
                </tr>
              </thead>
              <tbody className="divide-y divide-slate-100">
                <tr><td className="py-1.5 pr-3">Donahoe et al. (1993)</td><td className="py-1.5">{t.help.studyDonahoe1993}</td></tr>
                <tr><td className="py-1.5 pr-3">Burgos (1997)</td><td className="py-1.5">{t.help.studyBurgos1997}</td></tr>
                <tr><td className="py-1.5 pr-3">Donahoe & Burgos (1999)</td><td className="py-1.5">{t.help.studyDonahoeBurgos1999}</td></tr>
                <tr><td className="py-1.5 pr-3">Burgos (2000)</td><td className="py-1.5">{t.help.studyBurgos2000}</td></tr>
                <tr><td className="py-1.5 pr-3">Donahoe & Burgos (2000)</td><td className="py-1.5">{t.help.studyDonahoeBurgos2000}</td></tr>
                <tr><td className="py-1.5 pr-3">Burgos (2003)</td><td className="py-1.5">{t.help.studyBurgos2003}</td></tr>
                <tr><td className="py-1.5 pr-3">Burgos (2005)</td><td className="py-1.5">{t.help.studyBurgos2005}</td></tr>
                <tr><td className="py-1.5 pr-3">Burgos & Murillo-Rodr&iacute;guez (2007)</td><td className="py-1.5">{t.help.studyBurgosMurillo2007}</td></tr>
                <tr><td className="py-1.5 pr-3">Burgos (2007)</td><td className="py-1.5">{t.help.studyBurgos2007}</td></tr>
                <tr><td className="py-1.5 pr-3">Burgos et al. (2008)</td><td className="py-1.5">{t.help.studyBurgos2008}</td></tr>
                <tr><td className="py-1.5 pr-3">S&aacute;nchez et al. (2010)</td><td className="py-1.5">{t.help.studySanchez2010}</td></tr>
                <tr><td className="py-1.5 pr-3">Burns et al. (2011)</td><td className="py-1.5">{t.help.studyBurns2011}</td></tr>
                <tr><td className="py-1.5 pr-3">Calvin & McDowell (2015)</td><td className="py-1.5">{t.help.studyCalvin2015}</td></tr>
                <tr><td className="py-1.5 pr-3">Burgos & Donahoe (2016)</td><td className="py-1.5">{t.help.studyBurgosDonahoe2016}</td></tr>
                <tr><td className="py-1.5 pr-3">Alcal&aacute; (2017)</td><td className="py-1.5">{t.help.studyAlcala2017}</td></tr>
                <tr><td className="py-1.5 pr-3">Burgos (2019)</td><td className="py-1.5">{t.help.studyBurgos2019}</td></tr>
                <tr><td className="py-1.5 pr-3">Castiello et al. (2020)</td><td className="py-1.5">{t.help.studyCastiello2020}</td></tr>
                <tr><td className="py-1.5 pr-3">Burgos & Galeazzi (2021)</td><td className="py-1.5">{t.help.studyBurgosGaleazzi2021}</td></tr>
                <tr><td className="py-1.5 pr-3">Ojeda-Aguilar et al. (2023)</td><td className="py-1.5">{t.help.studyOjeda2023}</td></tr>
                <tr><td className="py-1.5 pr-3">Aguayo-Mendoza et al. (2024)</td><td className="py-1.5">{t.help.studyAguayo2024}</td></tr>
                <tr><td className="py-1.5 pr-3">Casta&ntilde;eda et al. (2025)</td><td className="py-1.5">{t.help.studyCastaneda2025}</td></tr>
              </tbody>
            </table>
          </div>
        </div>
      </Card>

      {/* Technical Support */}
      <Card className="p-6">
        <div className="flex items-center gap-2 mb-4">
          <BookOpen size={18} className="text-cyan-600" />
          <h2 className="text-lg font-bold text-slate-800">{t.help.supportTitle}</h2>
        </div>
        <div className="space-y-3 text-sm text-slate-600 leading-relaxed">
          <p>{t.help.supportIntro}</p>
          <ul className="ml-4 space-y-1">
            <li><span className="font-semibold text-slate-700">Miguel &Aacute;ngel Aguayo Mendoza</span></li>
            <li>{t.help.contactEmail}{' '}
              <a href="mailto:miguel.aguayo@academicos.udg.mx" className="text-cyan-600 hover:text-cyan-700 underline">
                miguel.aguayo@academicos.udg.mx
              </a>
            </li>
            <li>{t.help.contactUniversity}</li>
          </ul>
          <p>
            {t.help.labInfo}{' '}
            <a href="https://ceic.cucba.udg.mx/investigacion/laboratorios/investigacion-experimental" target="_blank" rel="noopener noreferrer"
              className="inline-flex items-center gap-1 text-cyan-600 hover:text-cyan-700 underline">
              <ExternalLink size={10} /> {t.help.website}
            </a>
          </p>
          <p className="text-xs text-slate-400 mt-2">
            {t.app.credits} &mdash; {t.app.lastModified}
          </p>
        </div>
      </Card>
    </div>
    </PageTransition>
  );
}
