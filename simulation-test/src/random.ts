// ══════════════════════════════════════════════════════════════════════════════
// DDM-UI — Random Number Generators
// Equivalent to R's rnorm() and rbeta() using IEEE 754 double precision
// ══════════════════════════════════════════════════════════════════════════════

/**
 * Seedable pseudo-random number generator (xoshiro128**)
 * Provides reproducible results when seed is set.
 * When no seed is provided, uses Math.random() as fallback.
 */
export class RNG {
  private state: Uint32Array | null = null;

  constructor(seed?: number) {
    if (seed !== undefined) {
      this.state = new Uint32Array(4);
      // SplitMix32 to initialize state from a single seed
      let s = seed >>> 0;
      for (let i = 0; i < 4; i++) {
        s = (s + 0x9e3779b9) >>> 0;
        let z = s;
        z = (z ^ (z >>> 16)) >>> 0;
        z = Math.imul(z, 0x85ebca6b) >>> 0;
        z = (z ^ (z >>> 13)) >>> 0;
        z = Math.imul(z, 0xc2b2ae35) >>> 0;
        z = (z ^ (z >>> 16)) >>> 0;
        this.state[i] = z;
      }
    }
  }

  /** Returns a uniform random number in [0, 1) */
  random(): number {
    if (this.state === null) return Math.random();

    const s = this.state;
    const result = Math.imul(s[1] * 5, 7) >>> 0;
    const t = (s[1] << 9) >>> 0;

    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];
    s[2] ^= t;
    s[3] = ((s[3] << 11) | (s[3] >>> 21)) >>> 0;

    return (result >>> 0) / 4294967296;
  }

  /**
   * Box-Muller transform: generates a standard normal variate.
   * Equivalent to R's rnorm(1, 0, 1)
   */
  rnorm(mean = 0, sd = 1): number {
    let u1 = this.random();
    let u2 = this.random();
    // Avoid log(0)
    while (u1 === 0) u1 = this.random();
    const z = Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
    return mean + sd * z;
  }

  /**
   * Beta distribution variate using Jöhnk's algorithm for small params
   * and the standard gamma-ratio method for larger params.
   * Equivalent to R's rbeta(1, alpha, beta)
   */
  rbeta(alpha: number, beta: number): number {
    const ga = this.rgamma(alpha);
    const gb = this.rgamma(beta);
    return ga / (ga + gb);
  }

  /**
   * Gamma distribution variate using Marsaglia & Tsang's method.
   * Equivalent to R's rgamma(1, shape, 1)
   */
  private rgamma(shape: number): number {
    if (shape < 1) {
      // For shape < 1, use the relation: Gamma(a) = Gamma(a+1) * U^(1/a)
      const u = this.random();
      return this.rgamma(shape + 1) * Math.pow(u, 1 / shape);
    }

    const d = shape - 1 / 3;
    const c = 1 / Math.sqrt(9 * d);

    while (true) {
      let x: number;
      let v: number;
      do {
        x = this.rnorm();
        v = 1 + c * x;
      } while (v <= 0);

      v = v * v * v;
      const u = this.random();

      if (u < 1 - 0.0331 * (x * x) * (x * x)) return d * v;
      if (Math.log(u) < 0.5 * x * x + d * (1 - v + Math.log(v))) return d * v;
    }
  }

  /**
   * Uniform random integer in [min, max] (inclusive).
   * Equivalent to R's sample(min:max, 1)
   */
  randInt(min: number, max: number): number {
    return min + Math.floor(this.random() * (max - min + 1));
  }

  /**
   * Fisher-Yates shuffle (in-place).
   * Equivalent to R's sample(x, length(x), replace=FALSE)
   */
  shuffle<T>(arr: T[]): T[] {
    for (let i = arr.length - 1; i > 0; i--) {
      const j = Math.floor(this.random() * (i + 1));
      [arr[i], arr[j]] = [arr[j], arr[i]];
    }
    return arr;
  }

  /** Returns a shuffled copy (does not modify original) */
  shuffled<T>(arr: T[]): T[] {
    const copy = [...arr];
    return this.shuffle(copy);
  }
}
