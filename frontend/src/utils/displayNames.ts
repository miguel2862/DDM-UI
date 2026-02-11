// Map digit characters to subscript equivalents
const subscripts: Record<string, string> = {
  '0': '\u2080', '1': '\u2081', '2': '\u2082', '3': '\u2083', '4': '\u2084',
  '5': '\u2085', '6': '\u2086', '7': '\u2087', '8': '\u2088', '9': '\u2089',
};

export function getDisplayName(name: string): string {
  // Special cases
  if (name === 'dVTA') return '\u03B4_VTA';
  if (name === 'dCA1') return '\u03B4_CA1';

  // Handle double dot notation: S..1 -> S''_1
  const doubleDotMatch = name.match(/^([A-Za-z]+)\.\.(\d+)$/);
  if (doubleDotMatch) {
    const [, base, num] = doubleDotMatch;
    const sub = num.split('').map(d => subscripts[d] || d).join('');
    return `${base}''\u200B${sub}`;
  }

  // Handle single dot notation: M.1 -> M'_1
  const singleDotMatch = name.match(/^([A-Za-z]+)\.(\d+)$/);
  if (singleDotMatch) {
    const [, base, num] = singleDotMatch;
    const sub = num.split('').map(d => subscripts[d] || d).join('');
    return `${base}'\u200B${sub}`;
  }

  // Handle plain name with trailing number: H1 -> H_1
  const plainNumMatch = name.match(/^([A-Za-z]+)(\d+)$/);
  if (plainNumMatch) {
    const [, base, num] = plainNumMatch;
    const sub = num.split('').map(d => subscripts[d] || d).join('');
    return `${base}${sub}`;
  }

  return name;
}
