export interface Phenomenon {
  id: string;
  name: string;
  description: string;
  category: string;
  available: boolean;
  phases: string;
  icon: string;
  color: string;
  npeCount: number;
  connectionCount: number;
}

export const phenomena: Phenomenon[] = [





  {
    id: 'extinction',
    name: 'Extinction',
    description: 'A conditioned response decreases when the CS is repeatedly presented without the US. The network learns to suppress responses that no longer predict reinforcement.',
    category: 'Basic',
    available: true,
    phases: 'Training (100) \u2192 Extinction (100)',
    icon: 'Flame',
    color: 'cyan',
    npeCount: 7,
    connectionCount: 6,
  },
  {
    id: 'alcala_2017',
    name: 'Autoshaped Impulsivity (Alcalá, 2017)',
    description: 'Smaller-Sooner (SS) vs Larger-Later (LL) autoshaped impulsive choice with ITI. Context stimulus modulates responding to immediate and delayed reinforcers.',
    category: 'Choice',
    available: true,
    phases: 'SS/LL (100 each, ITI) \u2192 Test (25)',
    icon: 'Scale',
    color: 'amber',
    npeCount: 13,
    connectionCount: 13,
  },
  {
    id: 'burgos_donahoe_blocking',
    name: 'Blocking (Burgos & Donahoe, 2016)',
    description: 'Kamin blocking: prior A+ training prevents learning about X when AX+ is presented. Demonstrates that redundant predictors fail to acquire associative strength.',
    category: 'Compound',
    available: true,
    phases: 'A+ (100) \u2192 AX+ (100) \u2192 Test X (25)',
    icon: 'ShieldOff',
    color: 'violet',
    npeCount: 13,
    connectionCount: 17,
  },
  {
    id: 'burgos_donahoe_successive',
    name: 'Successive (Burgos & Donahoe, 2016)',
    description: 'Successive training: A+ then X+ are trained separately. Both acquire associative strength independently. Baseline comparison for blocking.',
    category: 'Compound',
    available: true,
    phases: 'A+ (100) \u2192 X+ (100) \u2192 Test A (25)',
    icon: 'ArrowRightLeft',
    color: 'emerald',
    npeCount: 13,
    connectionCount: 17,
  },
  {
    id: 'acquisition',
    name: 'Acquisition',
    description: 'The most fundamental learning phenomenon: a neutral CS paired with a biologically significant US gradually acquires the ability to elicit a conditioned response.',
    category: 'Basic',
    available: true,
    phases: 'Training CS+US (100)',
    icon: 'Zap',
    color: 'cyan',
    npeCount: 7,
    connectionCount: 6,
  },
  {
    id: 'latent_inhibition',
    name: 'Latent Inhibition',
    description: 'Pre-exposure to a CS without consequence retards later conditioning when that CS is paired with a US. The prior non-reinforced experience makes the CS less associable.',
    category: 'Basic',
    available: true,
    phases: 'Pre-exposure CS (100) \u2192 Training CS+US (100)',
    icon: 'EyeOff',
    color: 'teal',
    npeCount: 7,
    connectionCount: 6,
  },
];
