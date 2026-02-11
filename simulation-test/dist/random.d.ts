/**
 * Seedable pseudo-random number generator (xoshiro128**)
 * Provides reproducible results when seed is set.
 * When no seed is provided, uses Math.random() as fallback.
 */
export declare class RNG {
    private state;
    constructor(seed?: number);
    /** Returns a uniform random number in [0, 1) */
    random(): number;
    /**
     * Box-Muller transform: generates a standard normal variate.
     * Equivalent to R's rnorm(1, 0, 1)
     */
    rnorm(mean?: number, sd?: number): number;
    /**
     * Beta distribution variate using Jöhnk's algorithm for small params
     * and the standard gamma-ratio method for larger params.
     * Equivalent to R's rbeta(1, alpha, beta)
     */
    rbeta(alpha: number, beta: number): number;
    /**
     * Gamma distribution variate using Marsaglia & Tsang's method.
     * Equivalent to R's rgamma(1, shape, 1)
     */
    private rgamma;
    /**
     * Uniform random integer in [min, max] (inclusive).
     * Equivalent to R's sample(min:max, 1)
     */
    randInt(min: number, max: number): number;
    /**
     * Fisher-Yates shuffle (in-place).
     * Equivalent to R's sample(x, length(x), replace=FALSE)
     */
    shuffle<T>(arr: T[]): T[];
    /** Returns a shuffled copy (does not modify original) */
    shuffled<T>(arr: T[]): T[];
}
