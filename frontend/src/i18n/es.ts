// Spanish translations for DDM-UI
import type { Translations } from './en';

export const es: Translations = {
  // ── App-wide ──
  app: {
    title: 'DDM-UI',
    subtitle: 'Discrepancia Temporal Difusa',
    subtitleFormer: 'anteriormente conocido como Modelo de Discrepancia Difusa (DiffDiscM)',
    tagline: 'Basado en Donahoe, Burgos y Palmer (1993)',
    taglineLong: 'Una interpretación conexionista del principio unificado de reforzamiento',
    credits: 'Interfaz diseñada por Miguel Ángel Aguayo Mendoza',
    university: 'Universidad de Guadalajara',
    lastModified: 'Última modificación: Septiembre 2026',
    version: 'v3.2',
  },

  // ── Navigation / Sidebar ──
  nav: {
    dashboard: 'Inicio',
    network: 'Red',
    trials: 'Ensayos',
    simulate: 'Simular',
    results: 'Resultados',
    paramSweep: 'Barrido',
    help: 'Ayuda',
    mode: 'Modo',
    beginner: 'Principiante',
    expert: 'Experto',
  },

  // ── Dashboard ──
  dashboard: {
    heroTitle1: 'Modelo de',
    heroTitleOld: 'Discrepancia Difusa',
    heroTitleNew: 'Discrepancia Temporal Difusa',
    heroDescription: 'Construya redes, defina protocolos temporales y simule el aprendizaje con el DDM original (Donahoe, Burgos y Palmer, 1993). Exporte redes para publicación e inspeccione las ecuaciones reales timestep por timestep.',
    quickStart: 'Inicio rápido: Extinción',
    buildScratch: 'Construir desde cero',

    npes: 'ENPs',
    connections: 'Conexiones',
    trialTypes: 'Tipos de ensayo',
    phases: 'Fases',
    architecture: 'Arquitectura',
    trialsStep: 'Ensayos',
    contingencies: 'Contingencias',
    simulateStep: 'Simular',
    resultsStep: 'Resultados',
  },

  // ── Gallery ──
  gallery: {
    title: 'Fenómenos de Condicionamiento',
    available: 'disponibles',
    ready: 'Listo',
    comingSoon: 'Próximamente',
    loadTemplate: 'Cargar plantilla',
    loading: 'Cargando...',
    templateNotAvailable: 'Plantilla aún no disponible',
    npesLabel: 'ENPs',
    connectionsLabel: 'conexiones',
  },

  // ── Phenomena (names + descriptions) ──
  phenomena: {





    extinction: {
      name: 'Extinción',
      description: 'La respuesta condicionada disminuye cuando el EC se presenta repetidamente sin el EI. La red aprende a suprimir respuestas que ya no predicen reforzamiento.',
      category: 'Básico',
      phases: 'Entrenamiento (100) → Extinción (100)',
    },
    alcala_2017: {
      name: 'Impulsividad Automoldeada (Alcalá, 2017)',
      description: 'Elección impulsiva automoldeada Menor-Pronto (SS) vs Mayor-Después (LL) con IEE. El estímulo de contexto modula la respuesta a reforzadores inmediatos y demorados.',
      category: 'Elección',
      phases: 'SS/LL (100 cada una, IEE) → Prueba (25)',
    },
    burgos_donahoe_blocking: {
      name: 'Bloqueo (Burgos & Donahoe, 2016)',
      description: 'Bloqueo de Kamin: el entrenamiento previo con A+ previene el aprendizaje sobre X cuando se presenta AX+. Demuestra que los predictores redundantes no adquieren fuerza asociativa.',
      category: 'Compuesto',
      phases: 'A+ (100) → AX+ (100) → Prueba X (25)',
    },
    burgos_donahoe_successive: {
      name: 'Sucesivo (Burgos & Donahoe, 2016)',
      description: 'Entrenamiento sucesivo: A+ y luego X+ se entrenan por separado. Ambos adquieren fuerza asociativa independientemente. Comparación de línea base para bloqueo.',
      category: 'Compuesto',
      phases: 'A+ (100) → X+ (100) → Prueba A (25)',
    },
    acquisition: {
      name: 'Adquisición',
      description: 'El fenómeno de aprendizaje más fundamental: un EC neutro pareado con un EI biológicamente significativo adquiere gradualmente la capacidad de elicitar una respuesta condicionada.',
      category: 'Básico',
      phases: 'Entrenamiento EC+EI (100)',
    },
    latent_inhibition: {
      name: 'Inhibición Latente',
      description: 'La pre-exposición a un EC sin consecuencia retarda el condicionamiento posterior cuando ese EC se parea con un EI. En el DTD, las presentaciones previas sin reforzamiento debilitan los pesos de conexión mediante el mecanismo de discrepancia, retardando la adquisición subsecuente.',
      category: 'Básico',
      phases: 'Preexposición EC (100) → Entrenamiento EC+EI (100)',
    },
  },

  // ── Network Builder ──
  network: {
    pageTitle: 'Arquitectura de Red',
    pageSubtitle: 'Construya y edite la arquitectura de su red neural',
    tabUnits: 'Unidades',
    tabConnections: 'Conexiones',
    export: 'Exportar',
    import: 'Importar',
    addNPE: 'Agregar ENP',
    name: 'Nombre',
    namePlaceholder: 'ej., S1, M.1, H1...',
    type: 'Tipo',
    excitatory: 'Excitatorio',
    inhibitory: 'Inhibitorio',
    layer: 'Capa',
    freeParameters: 'Parámetros libres',
    activation: 'Activación',
    tauSum: 'τ (Sum)',
    kappaDecay: 'κ (Decaim.)',
    muMean: 'μ (Media)',
    sigmaDev: 'σ (Desv)',
    logistic: 'Logística',
    addUnit: 'Agregar unidad',
    addConnection: 'Agregar conexión',
    source: 'Origen',
    target: 'Destino',
    weight: 'Peso',
    alphaIncrement: 'α (Incremento)',
    betaDecrement: 'β (Decremento)',
    alphaprimeInhibInc: "α' (Inc. Inhib.)",
    betaprimeInhibDec: "β' (Dec. Inhib.)",
    organize: 'Organizar',
    organizeTitle: 'Organizar nodos automáticamente por capa',
    lock: 'Fijar',
    locked: 'Fijado ✓',
    lockTitle: 'Fijar disposición actual para usar en Resultados',
    exportPNG: 'Exportar PNG',
    select: 'Seleccionar...',
    duplicateName: 'Ya existe un NPE con este nombre.',
    disconnected: 'Sin conexiones',
    disconnectedWarning: 'Este NPE no tiene conexiones. No participará en la simulación. Agregue al menos una conexión desde o hacia esta unidad.',
    layoutTraditional: 'Tradicional',
    layoutModern: 'Moderno',
    layoutTransduction: 'Transducción',
    stimulus: 'estímulo',
    response: 'respuesta',
  },

  // ── Trial Designer ──
  trial: {










    pageTitle: 'Diseñador de Ensayos y Contingencias',
    pageSubtitle: 'Defina tipos de ensayo y configure fases experimentales',
    addTrial: 'Agregar ensayo',
    addITI: 'Agregar IEE',
    createTrial: 'Crear ensayo',
    trialName: 'Nombre del ensayo',
    trialNameTooltip: 'Un ensayo representa un evento de aprendizaje individual. Cada ensayo tiene una secuencia de pasos de tiempo donde diferentes estímulos pueden estar activos y el aprendizaje puede ocurrir.',
    trialNamePlaceholder: 'ej., Entrenamiento, Extinción...',
    timesteps: 'Pasos de tiempo',
    timestepsTooltip: 'Cada paso de tiempo dentro de un ensayo define qué estímulos están activos (1.0 = presente, 0.0 = ausente) y si se aplica la regla de aprendizaje.',
    generate: 'Generar',
    stimulusActivations: 'Activaciones de estímulo por paso de tiempo:',
    fill: 'Llenar',
    fillTitle: 'Llenar todos los valores de estímulo con valores predeterminados',
    clear: 'Limpiar',
    clearTitle: 'Establecer todos los valores de estímulo en 0',
    learnOn: 'Aprend. Sí',
    learnOnTitle: 'Activar regla de aprendizaje en todos los pasos de tiempo',
    learnOff: 'Aprend. No',
    learnOffTitle: 'Desactivar regla de aprendizaje en todos los pasos de tiempo',
    tHeader: 't',
    learnHeader: 'Aprend.',
    updateTrial: 'Actualizar ensayo',
    cancelEdit: 'Cancelar edición',
    configureITI: 'Configurar IEE',
    itiDescription: 'Un IEE (Intervalo Entre Ensayos) es un período de descanso de un paso de tiempo insertado entre ensayos. Todos los valores de estímulo son típicamente 0 (sin estimulación) y el aprendizaje generalmente está desactivado.',
    itiName: 'Nombre del IEE',
    itiNameTooltip: 'Nombre para el tipo de ensayo IEE. Aparecerá en la configuración de fase cuando active el IEE.',
    itiNamePlaceholder: 'ej., IEE',
    npeActivations: 'Activaciones de ENP (paso de tiempo único):',
    enableLearning: 'Activar regla de aprendizaje',
    definedTrials: 'Ensayos definidos',
    timestepLabel: 'pasos de tiempo',
    timestepITI: '1 paso de tiempo (IEE)',
    edit: 'Editar',
    delete: 'Eliminar',
    configurePhase: 'Configurar fase',
    phaseName: 'Nombre de la fase',
    phaseNameTooltip: 'Una fase es un bloque de entrenamiento con tipos de ensayo específicos. Por ejemplo, \'Entrenamiento\' con ensayos EC+EI seguido de \'Extinción\' con ensayos solo EC.',
    phaseNamePlaceholder: 'ej., entrenamiento, extinción...',
    presentationOrder: 'Orden de presentación',
    presentationOrderTooltip: 'Controla la secuencia de ensayos dentro de una fase. Aleatorio: todos los tipos de ensayo se agrupan y barajan aleatoriamente (ej., si tiene 100 SS y 100 LL, los 200 ensayos aparecen en orden aleatorio). En bloque: todos los ensayos de un tipo se presentan consecutivamente antes de que comience el siguiente tipo — comúnmente usado para fases de prueba. Alternado: los tipos de ensayo se intercalan estrictamente en orden (A, B, A, B, ...) — comúnmente usado para entrenamiento con múltiples condiciones.',
    random: 'Aleatorio',
    inBulk: 'En bloque',
    alternated: 'Alternado',
    trialTypesLabel: 'Tipos de ensayo',
    noTrialsDefined: 'No se han definido ensayos aún',
    trialCounts: 'Cantidad de ensayos',
    trialCountsTooltip: 'Número de veces que cada tipo de ensayo se presenta en esta fase. Use un guión para separar cantidades para múltiples tipos de ensayo (ej., \'100-50\' significa 100 del primer tipo y 50 del segundo).',
    trialCountsPlaceholder: 'ej., 100 o 100-50',
    resetActivations: 'Reiniciar activaciones',
    resetActivationsTooltip: 'Cuando está marcado, todas las activaciones de las unidades se reinician a cero al inicio de cada ensayo (HasITI = False). Desmarque para usar un ensayo IEE en su lugar — las activaciones se mantienen entre ensayos a través del período de descanso IEE (HasITI = True).',
    minITI: 'IEE Mín',
    maxITI: 'IEE Máx',
    itiTrial: 'Ensayo IEE',
    noITITrialDefined: 'No se ha definido un ensayo IEE aún. Use la pestaña "Agregar IEE" para crear uno primero.',
    addPhase: 'Agregar fase',
    updatePhase: 'Actualizar fase',
    selectITITrial: 'Seleccione un ensayo IEE para habilitar esta fase.',
    phasesLabel: 'Fases',
    trials: 'ensayos',
    moveUp: 'Subir',
    moveDown: 'Bajar',
    reset: 'Reinicio',
  },

  // ── Simulation ──
  sim: {






























    pageTitle: 'Simulación',
    pageSubtitle: 'Ejecute la red neural y observe el aprendizaje en tiempo real',
    networks: 'Redes',
    thresholdLabel: 'Umbral',
    status: 'Estado',
    discCriterion: 'Crit. Discr.',
    thresholdPreset: 'Umbral de activación (θ)',
    thresholdPresetTooltip: 'Determina el tipo de distribución y los parámetros del umbral de activación para todos los ENP. Gaussiano usa una distribución normal; Beta restringe los umbrales a [0, 1]. DDM-UI Predeterminado: θ ~ N(0.2, 0.15). Donahoe et al. (1993): θ ~ N(0, 1), los parámetros originales publicados con mayor variabilidad en el umbral. Beta DDM-UI: θ ~ Beta(0.2, 0.15), umbrales acotados entre 0 y 1.',
    preset_gaussian_ddmui: 'Gaussiano',
    preset_gaussian_donahoe1993: 'Gaussiano (1993)',
    preset_beta_ddmui: 'Beta',
    updateProcedure: 'Procedimiento de actualización',
    updateProcedureTooltip: 'Controla el orden en que los ENP y conexiones se procesan durante el paso de aprendizaje.',
    asyncRandom: 'Asíncr. Aleat.',
    asyncRandomTooltip: 'Asíncrono + Aleatorio: Cada ENP actualiza sus pesos de conexión inmediatamente antes de que el siguiente sea procesado (aprendizaje en línea). El orden de los ENP y conexiones se baraja aleatoriamente en cada paso de tiempo. Este es el predeterminado y coincide con el comportamiento original del DiffDiscM.',
    asyncSequential: 'Asíncr. Secuenc.',
    asyncSequentialTooltip: 'Asíncrono + Secuencial: Cada ENP actualiza sus pesos de conexión inmediatamente antes de que el siguiente sea procesado. Los ENP y conexiones siempre se procesan en el orden fijo en que fueron definidos en la arquitectura de la red.',
    syncRandom: 'Síncr. Aleat.',
    syncRandomTooltip: 'Síncrono + Aleatorio: Primero se calculan todas las activaciones de los ENP, luego todos los pesos se actualizan simultáneamente (aprendizaje por lotes). El orden de procesamiento se baraja aleatoriamente en cada paso de tiempo.',
    syncSequential: 'Síncr. Secuenc.',
    syncSequentialTooltip: 'Síncrono + Secuencial: Primero se calculan todas las activaciones de los ENP, luego todos los pesos se actualizan simultáneamente. Los ENP y conexiones se procesan en el orden fijo en que fueron definidos.',
    discrepancyCriterion: 'Criterio de discrepancia',
    discrepancyTooltip: 'El criterio de discrepancia (disc) selecciona la rama de aprendizaje: con d ≥ disc se aplica la regla de incremento; con d < disc, la de decremento. Predeterminado: 0.0015. Los experimentos guardados conservan el valor configurado.',
    runSimulation: 'Ejecutar simulación',
    completeBeforeRunning: 'Complete la arquitectura de red, ensayos y contingencias antes de ejecutar.',
    simulating: 'Simulando...',
    processing: 'Procesando',
    networksSuffix: 'red(es)',
    durationConnector: 'en',
    completed: 'completada(s)',
    simulationFailed: 'Simulación fallida',
    apiNotRunning: 'Asegúrese de que la API de R esté ejecutándose en el puerto 8000.',
    tryAgain: 'Intentar de nuevo',
    simulationComplete: 'Simulación completada',
    runAgain: 'Ejecutar de nuevo',
    viewResults: 'Ver resultados',
    saveExperiment: 'Guardar experimento',
    loadExperiment: 'Cargar experimento',
    saveExperimentTooltip: 'Guarda la arquitectura de red actual, definiciones de ensayos y contingencias como un archivo JSON. Puede recargarlo después para reproducir el mismo experimento.',
    loadExperimentTooltip: 'Carga un archivo de configuración de experimento previamente guardado (JSON). Esto reemplazará la red, ensayos y contingencias actuales.',
    experimentSaved: 'Experimento guardado exitosamente!',
    experimentLoaded: 'Experimento cargado exitosamente!',
    invalidFile: 'Archivo de experimento inválido. Seleccione un archivo JSON de experimento DDM-UI válido.',
    validationNoUSD: 'Advertencia: No se encontró conexión US→D. La señal de reforzamiento podría no funcionar correctamente.',
    validationDisconnected: 'Advertencia: Las siguientes unidades no tienen conexiones y no participarán en la simulación: {units}',
    errorFallback: 'No se pudo conectar con la API de R',
    errorTimeout: 'La simulación agotó el tiempo de espera. El modelo puede ser demasiado grande ({networks} redes × {phases} fases). Intente reducir el número de redes o ensayos.',
    errorNetwork: 'No se puede alcanzar el motor de simulación. Asegúrese de que la API de R esté ejecutándose en el puerto 8000.',
    errorInternal: 'El motor de simulación encontró un error interno. Revise la arquitectura de la red para configuraciones inválidas (p. ej., conexión US→D faltante, NPEs desconectados).',
    cancelSimulation: 'Cancelar',
    cancelling: 'Cancelando...',
    etaRemaining: '~{time} restante',
    etaCalculating: 'Calculando...',
  },

  // ── Results ──
  results: {
    pageTitle: 'Explorador de Resultados',
    timestepsLabel: 'pasos de tiempo',
    individual: 'Individual',
    general: 'General',
    networkLabel: 'Red',
    chartType: 'Tipo de gráfico',
    activations: 'Activaciones',
    weights: 'Pesos',
    aggregate: 'Por Fase',
    learningSignals: 'Señales de aprendizaje',
    timestepPerPhase: 'Paso de tiempo por fase',
    units: 'Unidades',
    connectionsLabel: 'Conexiones',
    phase: 'Fase',
    measure: 'Medida',
    mean: 'Media',
    median: 'Mediana',
    noResults: 'Sin resultados de simulación',
    noResultsDescription: 'Ejecute una simulación primero para ver resultados aquí. Vaya a la pestaña Simular y haga clic en Ejecutar.',
    unitActivationsMeasure: 'Activaciones de unidades',
    unitActivationsOverTrials: 'Activaciones de unidades por ensayo',
    connectionWeightsMeasure: 'Pesos de conexiones',
    connectionWeightsOverTrials: 'Pesos de conexiones por ensayo',
    learningSignalsMeasure: 'Señales de aprendizaje',
    learningSignalsOverTrials: 'Señales de aprendizaje por ensayo',
    aggregateActivations: 'Activaciones agregadas',
    allNetworks: 'Todas las redes',
    thisNetwork: 'Esta red',
    allNetworksExport: 'Todas las redes',
    selected: 'Selección',
    downloadCurrent: 'Descargar red actual',
    downloadAll: 'Descargar todas las redes combinadas',
    downloadSelected: 'Descargar solo columnas seleccionadas',
    networksLegend: 'Redes:',
    net: 'Red',
    barsOverall: 'Barras = general',
    dotsPerNetwork: 'puntos = por red',
    dVTA: 'dVTA (Dopaminérgico)',
    dH: 'dH (Hipocampal)',
    // APA figure titles (for PNG export)
    figActivationsIndividual: 'Activaciones de Unidades por Ensayo',
    figActivationsGeneral: 'Activaciones de Unidades ({measure}) — Todas las Redes',
    figWeightsIndividual: 'Pesos de Conexiones por Ensayo',
    figWeightsGeneral: 'Pesos de Conexiones ({measure}) — Todas las Redes',
    figSignalsIndividual: 'Señales de Aprendizaje por Ensayo',
    figSignalsGeneral: 'Señales de Aprendizaje ({measure}) — Todas las Redes',
    figAggregate: 'Activación de Unidades ({measure}) Durante {phase}',
    figNote: '{networkInfo}. Pasos de tiempo: {timestepInfo}. Generado por DDM-UI v3.0.',
    figNetworkIndividual: 'Red {current} de {total}',
    figNetworkGeneral: '{total} redes ({measure})',
    // Extras
    exportPNG: 'Exportar gráfico como imagen PNG',
    selectAll: 'Seleccionar todo',
    deselectAll: 'Deseleccionar todo',
  },

  // ── Parameter Sweep ──
  sweep: {
    pageTitle: 'Barrido de Parámetros',
    pageSubtitle: 'Varíe sistemáticamente un parámetro y mida su efecto en la salida de la red',
    howItWorks: 'Cómo funciona el Barrido de Parámetros:',
    howDescription1: 'Seleccione una <strong>conexión o unidad</strong> y uno de sus parámetros (ej., peso, alfa). El barrido ejecutará la simulación completa múltiples veces, cada vez con un valor diferente de ese parámetro dentro del rango que especifique.',
    howDescription2: 'La <strong>Unidad de Salida</strong> es la unidad cuya activación se mide al final de cada simulación. Elija la unidad que desea observar (típicamente D para dopaminérgica o una unidad motora como M.1). El efecto puede ser <em>indirecto</em> — cambiar el peso de S1→H1 afecta a D a través de la cadena de conexiones de la red.',
    sweepConfig: 'Configuración del barrido',
    targetType: 'Tipo de objetivo',
    targetTypeTooltip: 'Elija si desea variar un parámetro de Conexión (ej., peso entre dos unidades) o un parámetro de Unidad (ej., media del umbral μ de un ENP).',
    connection: 'Conexión',
    unitNPE: 'Unidad (ENP)',
    element: 'Elemento',
    elementTooltip: 'La conexión o unidad específica cuyo parámetro será variado. Para conexiones, el formato es Origen-Destino (ej., S1-H1 significa la conexión de S1 a H1).',
    parameter: 'Parámetro',
    parameterTooltip: 'El parámetro a variar. Para conexiones: peso (fuerza inicial), alfa/beta (tasas de aprendizaje para incremento/decremento). Para unidades: mu (media del umbral), sigma (desviación del umbral), tau (sumación temporal), kappa (decaimiento de activación).',
    min: 'Mín',
    minTooltip: 'Valor mínimo del rango de parámetro a probar.',
    max: 'Máx',
    maxTooltip: 'Valor máximo del rango de parámetro a probar.',
    steps: 'Pasos',
    stepsTooltip: 'Número de valores equiespaciados a probar entre Mín y Máx. Más pasos = mayor resolución pero mayor tiempo de ejecución.',
    networksPerStep: 'Redes/paso',
    networksPerStepTooltip: 'Cuántas redes estocásticas promediar por valor de parámetro. Más redes = resultados más suaves pero más lento.',
    outputUnit: 'Unidad de salida',
    outputUnitTooltip: 'La unidad cuya activación media (en la última fase) se grafica. Elija la unidad que le interesa medir — típicamente la salida motora (M.1) o la unidad dopaminérgica (D). El efecto del parámetro barrido se propaga a través de la red hasta esta unidad.',
    runSweep: 'Ejecutar barrido',
    running: 'Ejecutando',
    sensitivityAnalysis: 'Análisis de Sensibilidad',
    points: 'puntos',
    meanLabel: 'Media',
    medianLabel: 'Mediana',
    emptyState: 'Seleccione un elemento, parámetro, rango y unidad de salida. El barrido ejecuta la simulación completa en cada paso, variando solo el parámetro seleccionado.',
    runningSweep: 'Ejecutando barrido de parámetros...',
    activationLabel: 'activación',
    select: 'Seleccionar...',
  },

  // ── Help ──
  help: {
    pageTitle: 'Ayuda e Información',
    pageSubtitle: 'Conozca el DDM y cómo usar esta aplicación',
    aboutTitle: 'Acerca del Modelo',
    aboutP1: 'El Modelo de Discrepancia Temporal Difusa (DTD), anteriormente conocido como Modelo de Discrepancia Difusa (DiffDiscM), fue desarrollado originalmente por Donahoe, Burgos y Palmer (1993) como una interpretación conexionista del principio unificado de reforzamiento. Proporciona un marco computacional que da cuenta tanto del condicionamiento operante como del pavloviano dentro de una arquitectura de red neural única.',
    aboutP2: 'El modelo simula cómo los elementos de procesamiento neural (ENP) interactúan a través de conexiones ponderadas para producir comportamiento aprendido. El reforzamiento ocurre cuando una discrepancia entre los resultados esperados y reales activa una unidad dopaminérgica, que a su vez modula difusamente los pesos de conexión en toda la red.',
    aboutP3: 'DDM-UI fue desarrollado por Aguayo-Mendoza y Dos Santos (2025) como una interfaz de usuario para facilitar el uso del modelo DTD en investigación conductual. La interfaz original SelNet1© se utilizó en más de 20 estudios publicados que abarcan fenómenos como adquisición, extinción, bloqueo, ensombrecimiento, inhibición latente, automoldeamiento y condicionamiento contextual.',








    networkTitle: 'Arquitectura de Red',
    networkIntro: 'La red del DTD está compuesta por siete tipos de capas, cada una representando un rol funcional distinto:',
    layerUS: 'EI (Estímulo Incondicionado)',
    layerUSDesc: 'Representa la entrada del estímulo incondicionado a la red. Se activa cuando se presenta un evento biológicamente significativo.',
    layerPS: 'SensorialPrimario (S)',
    layerPSDesc: 'Corresponde a la corteza sensorial primaria. Recibe entrada sensorial directa de estímulos condicionados (ej., tonos, luces).',
    layerAS: 'SensorialAsociativo (S″)',
    layerASDesc: 'Representa áreas de asociación sensorial de orden superior. Integra información sensorial y forma asociaciones entre estímulos.',
    layerH: 'Hipocampal (H)',
    layerHDesc: 'Modela la función hipocampal, involucrada en el procesamiento contextual y configural de estímulos.',
    layerAM: 'MotorAsociativo (M″)',
    layerAMDesc: 'Representa áreas de asociación premotora. Conecta el procesamiento sensorial con la planificación de la salida motora.',
    layerPM: 'MotorPrimario (M′)',
    layerPMDesc: 'Corresponde a la corteza motora primaria. Genera respuestas conductuales cuando se activa suficientemente.',
    layerD: 'Dopaminérgico (D)',
    layerDDesc: 'Modela el sistema dopaminérgico (ej., ATV). Genera una señal de reforzamiento basada en la discrepancia entre los resultados esperados y reales, que modula difusamente los pesos sinápticos.',
    namingTitle: 'Convención de Nomenclatura',
    namingIntro: 'Los nombres de ENP siguen una notación compacta que codifica su tipo de capa:',
    namingSdd: 'Representa S\u2033 (S doble prima), una unidad sensorial asociativa. Los dos puntos indican doble prima.',
    namingMdd: 'Representa M\u2033 (M doble prima), una unidad motora asociativa.',
    namingMd: 'Representa M\u2032 (M prima simple), una unidad motora primaria. Un solo punto indica prima simple.',
    namingH: 'Unidades hipocampales. En arquitecturas de dos vías, H1 se conecta a S\u20331 y H2 a S\u20332.',
    namingD: 'La unidad dopaminérgica.',
    namingUS: 'La unidad de estímulo incondicionado. La conexión EI\u2192D siempre tiene peso fijo = 1.0.',
    namingS: 'Unidades sensoriales primarias para diferentes estímulos.',
    freeParamsTitle: 'Parámetros Libres',
    freeParamsIntro: 'Cada ENP y conexión se rige por varios parámetros libres:',
    npeParams: 'Parámetros de ENP',
    connParams: 'Parámetros de Conexión',
    paramMu: 'Media de la distribución del umbral de activación. Predeterminado: 0.2.',
    paramSigma: 'Desviación estándar de la distribución del umbral de activación. Predeterminado: 0.15.',
    paramTau: 'Controla cómo se acumula la activación en pasos de tiempo consecutivos. Predeterminado: 0.1.',
    paramKappa: 'Determina qué tan rápido disminuye la activación de un ENP cuando cesa la entrada. Predeterminado: 0.1.',
    paramLogistic: 'Pendiente de la función de activación logística. Predeterminado: 0.1.',
    paramWeight: 'Peso sináptico de una conexión. Predeterminado inicial: 0.1. EI\u2192D es fijo en 1.0.',
    paramAlpha: 'Tasa de aprendizaje para incrementos de peso excitatorio. Predeterminado: 0.5.',
    paramBeta: 'Tasa de aprendizaje para decrementos de peso excitatorio. Predeterminado: 0.12.',
    paramAlphaPrime: 'Tasa de aprendizaje para incrementos de peso inhibitorio. Predeterminado: 0.5.',
    paramBetaPrime: 'Tasa de aprendizaje para decrementos de peso inhibitorio. Predeterminado: 0.12.',
    howToUseTitle: 'Cómo Usar',
    step1Title: 'Cargue una plantilla o construya desde cero.',
    step1Desc: 'Desde el Inicio, seleccione un fenómeno preconstruido (ej., extinción, bloqueo) o navegue a la página de Red para definir sus propios ENP y conexiones.',
    step2Title: 'Configure ensayos.',
    step2Desc: 'En la página de Ensayos, defina tipos de ensayo y sus secuencias de pasos de tiempo. Cada paso de tiempo especifica qué ENP reciben entrada externa y si la regla de aprendizaje está activa.',
    step3Title: 'Establezca contingencias y fases.',
    step3Desc: 'Organice los tipos de ensayo en fases con números específicos de repeticiones. Puede agregar múltiples fases para modelar adquisición, extinción, períodos de descanso y condiciones de prueba.',
    step4Title: 'Ejecute la simulación.',
    step4Desc: 'Vaya a Simular, elija el número de redes y ejecute la simulación. El backend de R utiliza el DDM original. Opcionalmente, active el inspector de ecuaciones antes de ejecutar para registrar la primera red.',
    step5Title: 'Analice resultados.',
    step5Desc: 'En la página de Resultados, explore trayectorias de activación y pesos a través de ensayos y fases. Alterne entre vistas Individual (por red) y General (todas las redes). Use el gráfico agregado para comparar activaciones de unidades por fase.',
    referencesTitle: 'Referencias',
    studiesTitle: 'Estudios que Usan el Modelo DTD',
    studiesIntro: 'La siguiente es una selección de estudios publicados que han empleado la interfaz SelNet1© y/o el modelo DTD (anteriormente DiffDiscM) para investigación conductual computacional:',
    studyCol: 'Estudio',
    phenomenonCol: 'Fenómeno',
    // Studies table - phenomenon descriptions
    studyDonahoe1993: 'Adquisición, extinción, readquisición',
    studyBurgos1997: 'Redes neuronales artificiales en ambientes pavlovianos',
    studyDonahoeBurgos1999: 'Temporalidad sin temporizador',
    studyBurgos2000: 'Conducta supersticiosa',
    studyDonahoeBurgos2000: 'Re-evaluación del reforzamiento',
    studyBurgos2003: 'Inhibición latente',
    studyBurgos2005: 'Efectos de la razón C/T',
    studyBurgosMurillo2007: 'Especificidad contextual y renovación',
    studyBurgos2007: 'Automoldeamiento y automantenimiento',
    studyBurgos2008: 'Condicionamiento simultáneo',
    studySanchez2010: 'Condicionamiento de segundo orden',
    studyBurns2011: 'Condicionamiento pavloviano (membrana nictitante)',
    studyCalvin2015: 'Efecto de bloqueo',
    studyBurgosDonahoe2016: 'Mecanismos de bloqueo y ensombrecimiento',
    studyAlcala2017: 'Impulsividad pavloviana',
    studyBurgos2019: 'Visión ciega pavloviana y condicionamiento',
    studyCastiello2020: 'Contraste conductual pavloviano',
    studyBurgosGaleazzi2021: 'Rol hipocampal en condicionamiento pavloviano',
    studyOjeda2023: 'Ensombrecimiento post-entrenamiento',
    studyAguayo2024: 'Impulsividad automoldeada',
    studyCastaneda2025: 'Vías convergentes vs independientes en COS y bloqueo',
    // Contact
    contactEmail: 'Correo:',
    contactUniversity: 'Universidad de Guadalajara',
    supportTitle: 'Soporte Técnico',
    supportIntro: 'Si tiene algún problema o preguntas adicionales, por favor contacte a:',
    labInfo: 'Para más información, visite el laboratorio de investigación experimental y teórica en Aprendizaje, Condicionamiento y Conducta Adaptativa:',
    website: 'Sitio web',
  },

  // ── Session Restore ──
  session: {
    restoreTitle: 'Sesión anterior encontrada',
    restoreDescription: '¿Desea continuar donde lo dejó?',
    restore: 'Restaurar sesión',
    newSession: 'Nueva sesión',
  },

  // ── Splash Screen ──
  splash: {
    versionHistory: 'Historial de versiones',
    loading: 'Iniciando motor de simulación...',
  },

  // ── Notifications ──
  toast: {
    autoSaved: 'Progreso guardado automáticamente',
    autoRestored: 'Sesión anterior restaurada',
    connectionLost: 'Conexión con el motor de simulación perdida',
    connectionRestored: 'Conexión restaurada',
    simulationError: 'Simulación fallida. Revise la arquitectura de red e intente de nuevo.',
    apiTimeout: 'Tiempo de espera agotado. El motor de simulación puede estar ocupado.',
    unexpectedError: 'Ocurrió un error inesperado',
    copied: 'Copiado al portapapeles',
    exportPNGError: 'Error al exportar gráfico como PNG. Intente de nuevo.',
    simulationCancelled: 'Simulación cancelada',
    navigationBlocked: 'La navegación está bloqueada mientras la simulación se ejecuta. Cancele o espere a que termine.',
  },

  // ── Confirmation dialogs ──
  confirm: {
    resetTitle: 'Reiniciar Configuración',
    resetMessage: 'Esto eliminará todos los ENP, conexiones, ensayos y fases. Esta acción no se puede deshacer.',
    resetConfirm: 'Reiniciar todo',
    deleteNPETitle: 'Eliminar ENP',
    deleteNPEMessage: 'Eliminar este ENP también eliminará todas sus conexiones. ¿Continuar?',
    deleteTrialTitle: 'Eliminar Ensayo',
    deleteTrialMessage: 'Este ensayo puede estar siendo usado en contingencias de fase. ¿Continuar?',
    deletePhaseTitle: 'Eliminar Fase',
    deletePhaseMessage: '¿Está seguro de que desea eliminar esta fase?',
    thresholdChangeTitle: 'Cambiar preset de umbral',
    thresholdChangeMessage: 'Cambiar el preset de umbral sobrescribirá los valores de μ (media) y σ (desviación estándar) de todos los ENP con los valores predeterminados del preset. Esta acción no se puede deshacer.',
    loadTemplateTitle: 'Cargar plantilla',
    loadTemplateMessage: 'Cargar una nueva plantilla reemplazará la arquitectura de red, ensayos y contingencias actuales. Los resultados de simulación existentes se perderán.',
    confirm: 'Confirmar',
    cancel: 'Cancelar',
  },

  // ── Error Boundary ──
  errorBoundary: {
    title: 'Algo salió mal',
    description: 'Ocurrió un error inesperado en la aplicación. Sus datos han sido guardados automáticamente.',
    reload: 'Recargar Aplicación',
    reset: 'Reiniciar y Recargar',
    details: 'Detalles del Error',
  },

  // ── Common / Shared ──
  common: {
    select: 'Seleccionar...',
    cancel: 'Cancelar',
    save: 'Guardar',
    edit: 'Editar',
    delete: 'Eliminar',
    close: 'Cerrar',
    language: 'Idioma',
  },
};
