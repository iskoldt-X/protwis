/**
 * DataTablesSVG.js
 * ----------------
 * Lightweight SVG export for plain HTML tables (Bootstrap-friendly).
 *
 * Highlights
 * - White (or auto) THEAD background.
 * - Zebra TBODY (grey/white) for consistent look across browsers.
 * - Thin inner grid, solid black outer frame.
 * - Optional special vertical divider (e.g. split between “Classes”/“Receptors”).
 * - Two header separator lines (after each THEAD row), drawn on top.
 * - Soft text wrapping with optional “prefer parentheses” split.
 * - Optional link handling: keep as clickable <a>, or strip to plain text.
 *
 * Usage
 *   DataTablesSVG.export('#gpcrTable', {
 *     filename: 'GPCRdb_Human_GPCR_Classes.svg',
 *     links: 'strip',                // 'strip' (default) | 'keep'
 *     linkColor: '#0b63c0',          // when links:'keep'
 *     linkUnderline: false,          // when links:'keep'
 *     headerBackground: 'white',     // 'white' (default) | 'auto'
 *     bodyFallback: '#f8f9fa',       // only used when zebraOff
 *     preferParenWrap: true,         // split "(...)" onto two lines when it fits
 *     zebraOn: true,                 // default true
 *     zebraColors: ['#f8f9fa','#fff'],
 *     headerSeparators: true,        // draw lines after each THEAD row
 *     specialDividerAt: 3            // x = left edge of col 4 (0-based index)
 *   });
 */
(function (global) {
  const NS    = 'http://www.w3.org/2000/svg';
  const XLINK = 'http://www.w3.org/1999/xlink';

  /* --------------------------- utilities --------------------------- */

  function triggerDownload(text, filename, type = 'image/svg+xml') {
    const blob = new Blob([text], { type });
    const url  = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    a.remove();
    URL.revokeObjectURL(url);
  }

  function stampedBase(base, ext) {
    const d = new Date(), pad = n => String(n).padStart(2,'0');
    const stamp = `${d.getFullYear()}-${pad(d.getMonth()+1)}-${pad(d.getDate())}_${pad(d.getHours())}-${pad(d.getMinutes())}-${pad(d.getSeconds())}`;
    return `${base}-${stamp}.${ext}`;
  }

  // Normalize any CSS color into rgba(...) so we can detect transparency
  function pxColor(css) {
    if (!css) return 'rgba(0,0,0,0)';
    if (/^rgba?\(/i.test(css)) return css;
    const tmp = document.createElement('div');
    tmp.style.background = css;
    document.body.appendChild(tmp);
    const c = getComputedStyle(tmp).backgroundColor || 'rgba(0,0,0,0)';
    tmp.remove();
    return c;
  }

  function isTransparent(c) {
    const m = /^rgba?\(([^)]+)\)/i.exec(c || '');
    if (!m) return false; // named colors -> treat as opaque
    const parts = m[1].split(',').map(s => parseFloat(s));
    return parts.length === 4 ? (parts[3] <= 0.001) : false;
  }

  // Simple SVG line helper (crisp with .5 coordinates)
  function line(x1, y1, x2, y2, stroke) {
    const p = document.createElementNS(NS, 'path');
    p.setAttribute('d', `M ${x1} ${y1} L ${x2} ${y2}`);
    p.setAttribute('stroke', stroke);
    p.setAttribute('fill', 'none');
    return p;
  }

  /* ----------------------- measuring + wrapping -------------------- */

  // Shared invisible SVG for measuring text width
  let MEASURE = document.getElementById('DTSVG_MEASURE_MIN');
  if (!MEASURE) {
    MEASURE = document.createElementNS(NS, 'svg');
    MEASURE.id = 'DTSVG_MEASURE_MIN';
    MEASURE.style.cssText = 'position:fixed;left:-10000px;top:0;opacity:0;pointer-events:none';
    document.body.appendChild(MEASURE);
  }

  function measureWidth(str, font) {
    const t = document.createElementNS(NS, 'text');
    t.textContent = str || '';
    t.setAttribute('font-family', font?.fontFamily || 'system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial');
    if (font?.size)   t.setAttribute('font-size', font.size);
    if (font?.weight) t.setAttribute('font-weight', font.weight);
    MEASURE.appendChild(t);
    const w = t.getComputedTextLength() || 0;
    t.remove();
    return w;
  }

  // Prefer splitting "(...)" into two lines when it fits
  function splitParenTwoLines(text) {
    const m = String(text || '').match(/^(.*?)(\s*\(([^)]+)\))\s*$/);
    if (!m) return null;
    return [m[1].trim(), `(${m[3].trim()})`];
  }

  // Greedy word wrap into lines that fit maxW (SVG units/pixels)
  function wrapToWidth(text, maxW, font, extras) {
    if (!text) return [''];
    const preferParen = !!(extras && extras.preferParenWrap);
    const par = splitParenTwoLines(text);
    if (par && preferParen) {
      const w1 = measureWidth(par[0], font);
      const w2 = measureWidth(par[1], font);
      if (w1 <= maxW && w2 <= maxW) return par;
      if (w1 <= maxW * 1.05 && w2 <= maxW * 1.05) return par; // small grace
    }
    const words = String(text).split(/\s+/);
    const lines = [];
    let cur = '';
    for (const w of words) {
      const next = cur ? `${cur} ${w}` : w;
      if (measureWidth(next, font) <= maxW || !cur) cur = next;
      else { lines.push(cur); cur = w; }
    }
    if (cur) lines.push(cur);
    return lines.length ? lines : [''];
  }

  /* --------------------------- SVG helpers ------------------------- */

  function rect(x, y, w, h, fill, stroke) {
    const r = document.createElementNS(NS, 'rect');
    r.setAttribute('x', x);
    r.setAttribute('y', y);
    r.setAttribute('width',  w);
    r.setAttribute('height', h);
    if (fill)   r.setAttribute('fill',   fill);
    if (stroke) r.setAttribute('stroke', stroke);
    return r;
  }

  function textNode(x, y, str, opts = {}) {
    const t = document.createElementNS(NS, 'text');
    t.setAttribute('x', x);
    t.setAttribute('y', y);
    t.setAttribute('text-anchor', opts.anchor || 'middle');
    t.setAttribute('dominant-baseline', opts.baseline || 'middle');
    t.setAttribute('font-family', opts.fontFamily || 'system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial');
    if (opts.size)      t.setAttribute('font-size', opts.size);
    if (opts.weight)    t.setAttribute('font-weight', opts.weight);
    if (opts.fill)      t.setAttribute('fill', opts.fill);
    if (opts.underline) t.setAttribute('text-decoration', 'underline');
    t.textContent = str;
    return t;
  }

  // Column widths from first tbody row (or last header row as fallback)
  function measureCols(table) {
    const tbody = table.querySelector('tbody');
    const thead = table.querySelector('thead');
    const probe = (tbody && tbody.rows[0]) || (thead && thead.rows[thead.rows.length - 1]);
    const out = [];
    if (!probe) return out;
    for (let c = 0; c < probe.cells.length; c++) {
      out.push(Math.round(probe.cells[c].getBoundingClientRect().width));
    }
    return out;
  }

  // Extract text (+link when keeping links)
  function extract(el, mode) {
    const text = (el.textContent || '').replace(/\s+/g, ' ').trim();
    if (mode !== 'keep') return { text };
    const a = el.querySelectorAll('a[href]');
    if (a.length === 1) return { text, href: a[0].getAttribute('href') || '' };
    return { text };
  }

  /* ---------------------------- builder ---------------------------- */

  function buildSVG(table, opts = {}) {
    const thead = table.querySelector('thead');
    const tbody = table.querySelector('tbody');
    if (!thead || !tbody) throw new Error('DataTablesSVG: expected <thead> and <tbody>.');

    // Column metrics + row heights
    const colW        = measureCols(table);
    const probe       = tbody.rows[0] || thead.rows[0];
    const rowH        = Math.round(probe.getBoundingClientRect().height);
    const headHeights = Array.from(thead.rows).map(tr => Math.round(tr.getBoundingClientRect().height));
    const headH       = headHeights.reduce((a, b) => a + b, 0);

    // Y positions for header separator lines (after each THEAD row)
    const headerLineYs = [];
    if (headHeights.length) {
      let acc = 0;
      for (let i = 0; i < headHeights.length; i++) { acc += headHeights[i]; headerLineYs.push(acc); }
    }

    // Canvas size
    const width  = colW.reduce((a, b) => a + b, 0) + 1;
    const height = headH + (tbody.rows.length * rowH) + 1;

    // Root SVG
    const svg = document.createElementNS(NS, 'svg');
    svg.setAttribute('xmlns', NS);
    svg.setAttribute('width',  width);
    svg.setAttribute('height', height);
    svg.setAttribute('viewBox', `0 0 ${width} ${height}`);

    // Cumulative Xs for fast cell placement
    const cumX = [0]; for (let i = 0; i < colW.length; i++) cumX.push(cumX[i] + colW[i]);

    /* ----- options & look ----- */

    const linkMode         = (opts.links || 'strip').toLowerCase(); // 'strip' | 'keep'
    const linkColor        = opts.linkColor || '#0b63c0';
    const linkUnderline    = !!opts.linkUnderline;

    const borderLight      = '#dddddd';  // inner grid
    const borderHard       = '#000000';  // outer frame + special dividers

    const headMode         = (opts.headerBackground || 'white'); // 'white' | 'auto'
    const bodyFallback     = opts.bodyFallback || '#f8f9fa';
    const preferParenWrap  = !!opts.preferParenWrap;

    const zebraOn          = opts.zebraOn ?? true;
    const zebraColors      = opts.zebraColors || ['#f8f9fa', '#ffffff'];

    const headerSeparators = opts.headerSeparators ?? true;

    // Draw a vertical black divider at the LEFT edge of column (index)
    // For your GPCR table: 3 → between “CLASSES” (3 cols) and “RECEPTORS” (rest)
    const specialDividerAt = (typeof opts.specialDividerAt === 'number') ? opts.specialDividerAt : 3;

    // Typography
    const fontHeader = { size: 14, weight: 600 };
    const fontBody   = { size: 14, weight: 500 };

    /* --------------------------- THEAD ------------------------------ */

    let y = 0;
    for (let r = 0; r < thead.rows.length; r++) {
      const tr = thead.rows[r];
      const h  = headHeights[r];
      let ci = 0;

      Array.from(tr.cells).forEach(th => {
        const span = th.colSpan || 1;
        const x1 = cumX[ci], x2 = cumX[ci + span], w = x2 - x1;

        // Header background: force white (or auto-detect)
        let bg;
        if (headMode === 'white') bg = '#ffffff';
        else {
          bg = pxColor(getComputedStyle(th).backgroundColor);
          if (isTransparent(bg)) bg = '#ffffff';
        }
        svg.appendChild(rect(x1, y, w, h, bg, borderLight));

        // Header text (wrap, centered)
        const { text, href } = extract(th, linkMode);
        if (text) {
          const maxW = w - 12; // padding
          const lines = wrapToWidth(text, maxW, fontHeader, { preferParenWrap });
          const totalH = lines.length * fontHeader.size;
          let yy = y + h/2 - totalH/2 + fontHeader.size/2;
          lines.forEach(lineStr => {
            const node = textNode(x1 + w/2, yy, lineStr, {
              size: fontHeader.size, weight: fontHeader.weight,
              fill: href ? linkColor : '#000', underline: href && linkUnderline
            });
            if (href) {
              const a = document.createElementNS(NS, 'a');
              a.setAttributeNS(XLINK, 'xlink:href', href);
              a.setAttribute('target', '_blank');
              a.appendChild(node);
              svg.appendChild(a);
            } else {
              svg.appendChild(node);
            }
            yy += fontHeader.size;
          });
        }
        ci += span;
      });

      y += h;
    }

    /* --------------------------- TBODY ------------------------------ */

    for (let r = 0; r < tbody.rows.length; r++) {
      const tr  = tbody.rows[r];
      const top = headH + r * rowH;

      for (let c = 0; c < tr.cells.length; c++) {
        const td = tr.cells[c];
        const x  = cumX[c], w = colW[c];

        // Background fill: zebra first (for consistency)
        let bg;
        if (zebraOn) {
          bg = zebraColors[r % 2];
        } else {
          // Legacy: try to read row/td background, then fallback
          const bgTr = pxColor(getComputedStyle(tr).backgroundColor);
          const bgTd = pxColor(getComputedStyle(td).backgroundColor);
          if (!isTransparent(bgTr))      bg = bgTr;
          else if (!isTransparent(bgTd)) bg = bgTd;
          else                            bg = bodyFallback;
        }
        svg.appendChild(rect(x, top, w, rowH, bg, borderLight));

        // Cell text
        const { text, href } = extract(td, linkMode);
        if (text) {
          const maxW  = w - 12;
          const lines = wrapToWidth(text, maxW, fontBody, { preferParenWrap });
          const totalH = lines.length * fontBody.size;
          let yy = top + rowH/2 - totalH/2 + fontBody.size/2;
          lines.forEach(lineStr => {
            const node = textNode(x + w/2, yy, lineStr, {
              size: fontBody.size, weight: fontBody.weight,
              fill: href ? linkColor : '#000', underline: href && linkUnderline
            });
            if (href) {
              const a = document.createElementNS(NS, 'a');
              a.setAttributeNS(XLINK, 'xlink:href', href);
              a.setAttribute('target', '_blank');
              a.appendChild(node);
              svg.appendChild(a);
            } else {
              svg.appendChild(node);
            }
            yy += fontBody.size;
          });
        }
      }
    }

    /* ------------ overlays drawn ON TOP of body cells --------------- */

    // Header separator lines (after each THEAD row)
    if (headerSeparators && headerLineYs.length) {
      headerLineYs.forEach(ypos => {
        svg.appendChild(line(0.5, ypos + 0.5, width - 0.5, ypos + 0.5, borderHard));
      });
    }

    // Special vertical divider (e.g., after 3rd col)
    if (Number.isFinite(specialDividerAt) && specialDividerAt >= 0 && specialDividerAt < cumX.length) {
      const x = cumX[specialDividerAt] + 0.5;
      svg.appendChild(line(x, 0.5, x, height - 0.5, borderHard));
    }

    // Outer border (black frame)
    const border = document.createElementNS(NS, 'path');
    border.setAttribute('d', `M 0.5 0.5 H ${width - 0.5} V ${height - 0.5} H 0.5 Z`);
    border.setAttribute('stroke', borderHard);
    border.setAttribute('fill', 'none');
    svg.appendChild(border);

    return svg;
  }

  /* --------------------------- public API -------------------------- */

  function exportSVG(selector, options = {}) {
    const el = (typeof selector === 'string') ? document.querySelector(selector) : selector;
    if (!el) return console.error('DataTablesSVG: selector not found:', selector);
    const svg = buildSVG(el, options);
    const xml = new XMLSerializer().serializeToString(svg);
    const base = (typeof selector === 'string' ? (selector.match(/#([\w-]+)/)?.[1] || 'table') : 'table');
    const filename = options.filename || stampedBase(base, 'svg');
    triggerDownload(xml, filename);
  }

  global.DataTablesSVG = { export: exportSVG };
})(window);
