import type { NPERow, ConnectionRow } from './types.js';
export interface TemplateData {
    id: string;
    name: string;
    npes: NPERow[];
    connections: ConnectionRow[];
    trials: Record<string, string[]>;
    contingencies: string[];
    hasITI: boolean[];
}
export declare function getExtinctionTemplate(): TemplateData;
export declare function getAcquisitionTemplate(): TemplateData;
export declare function getSpontaneousRecoveryTemplate(): TemplateData;
export declare function getLatentInhibitionTemplate(): TemplateData;
export declare function getBlockingTemplate(): TemplateData;
/** All templates by ID */
export declare const TEMPLATES: Record<string, () => TemplateData>;
