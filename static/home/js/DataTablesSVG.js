/**
 * DataTablesSVG.js
 * ----------------
 * Export the CURRENT DOM of an HTML table (#ccsTable) as a clean, editable SVG.
 * - Preserves header/row labels, numbers, cell background colors, borders
 * - Uses the computed cell sizes & font sizes from your CSS (so it matches the view)
 *
 * Usage:
 *   DataTablesSVG.export('#ccsTable','ccs-table')
 */

(function (global) {
  function makeFilename(base, ext) {
    const d = new Date(), pad = n => String(n).padStart(2, '0');
    const stamp = `${d.getFullYear()}-${pad(d.getMonth()+1)}-${pad(d.getDate())}_` +
                  `${pad(d.getHours())}-${pad(d.getMinutes())}-${pad(d.getSeconds())}`;
    return `${base || 'table'}-${stamp}.${ext}`;
  }

  function triggerDownload(text, filename, type='image/svg+xml') {
    const blob = new Blob([text], { type });
    const url  = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url; a.download = filename;
    document.body.appendChild(a); a.click(); a.remove();
    URL.revokeObjectURL(url);
  }

  // Extract two-line label from a <th> like "Class A<br><small>(Rhodopsin)</small>"
  function parseHeaderCell(th) {
    // Try to preserve the two-line structure you render
    const small = th.querySelector('small');
    const line1 = th.cloneNode(true);
    // Remove <small> from line1 clone
    Array.from(line1.querySelectorAll('small')).forEach(el => el.remove());
    const primary = line1.textContent.trim().replace(/\s+/g, ' ');
    const secondary = small ? small.textContent.trim().replace(/\s+/g, ' ') : '';
    return { primary, secondary };
  }

  function cssPx(el, prop) {
    const v = getComputedStyle(el).getPropertyValue(prop);
    const n = parseFloat(v);
    return isFinite(n) ? n : null;
  }

  function toRGB(color) {
    // Already rgb()/rgba()? return as is; otherwise let the browser normalize it
    if (/^rgba?\(/i.test(color)) return color;
    const tmp = document.createElement('div');
    tmp.style.color = color;
    document.body.appendChild(tmp);
    const rgb = getComputedStyle(tmp).color;
    document.body.removeChild(tmp);
    return rgb || 'rgb(255,255,255)';
  }

  function textYCentered(y, h, fontSizePx) {
    // SVG text is baseline-aligned; shift so text visually centers in the cell.
    // 0.35 is a good baseline correction for typical fonts.
    return y + (h + fontSizePx * (1 - 0.35)) / 2;
  }

  function exportTable(selector, baseName) {
    const table = (typeof selector === 'string') ? document.querySelector(selector) : selector;
    if (!table) {
      console.error('DataTablesSVG: table not found:', selector);
      return;
    }

    const thead = table.querySelector('thead');
    const tbody = table.querySelector('tbody');
    if (!thead || !tbody) {
      console.error('DataTablesSVG: expected <thead> and <tbody>.');
      return;
    }

    // Dimensions from the live table (matches your CSS variables)
    const firstRow = tbody.rows[0] || thead.rows[0];
    const firstDataCell = tbody.querySelector('td') || table.querySelector('td');
    const firstRowHeader = tbody.querySelector('th') || table.querySelector('th');

    if (!firstRow || !firstDataCell || !firstRowHeader) {
      console.error('DataTablesSVG: missing cells to infer sizes.');
      return;
    }

    const cellW = Math.round(firstDataCell.getBoundingClientRect().width);
    const cellH = Math.round(firstDataCell.getBoundingClientRect().height);
    const rowHeaderW = Math.round(firstRowHeader.getBoundingClientRect().width);
    const borderPx = 1; // your CSS uses 1px borders

    // Fonts (use computed values from the actual cells)
    const numFontSize = parseFloat(getComputedStyle(firstDataCell).fontSize) || 14;
    const headFontSize = parseFloat(getComputedStyle(thead.querySelector('th')).fontSize) || 16;
    // Secondary header is stored in <small>
    const smallEl = thead.querySelector('small') || table.querySelector('small');
    const head2FontSize = smallEl ? parseFloat(getComputedStyle(smallEl).fontSize) : Math.max(11, headFontSize - 4);
    const numFontWeight  = getComputedStyle(firstDataCell).fontWeight || 550;

    // Table shape
    const cols = thead.rows[0] ? thead.rows[0].cells.length : (firstRow.cells.length);
    const rows = tbody.rows.length + 1; // +1 for header row
    const width  = rowHeaderW + (cols - 1) * cellW + borderPx; // rightmost stroke space
    const height = (1 * cellH) + (tbody.rows.length) * cellH + borderPx;

    // Build SVG
    const svgNS = 'http://www.w3.org/2000/svg';
    const svg = document.createElementNS(svgNS, 'svg');
    svg.setAttribute('xmlns', svgNS);
    svg.setAttribute('width',  String(width));
    svg.setAttribute('height', String(height));
    svg.setAttribute('viewBox', `0 0 ${width} ${height}`);

    // Background
    const bg = document.createElementNS(svgNS, 'rect');
    bg.setAttribute('x', '0'); bg.setAttribute('y', '0');
    bg.setAttribute('width', String(width)); bg.setAttribute('height', String(height));
    bg.setAttribute('fill', '#ffffff');
    svg.appendChild(bg);

    // Helpers to add shapes/text
    function rect(x,y,w,h, fill, stroke) {
      const r = document.createElementNS(svgNS, 'rect');
      r.setAttribute('x', x); r.setAttribute('y', y);
      r.setAttribute('width', w); r.setAttribute('height', h);
      if (fill)   r.setAttribute('fill', fill);
      if (stroke) r.setAttribute('stroke', stroke);
      return r;
    }
    function text(x,y, content, opts={}) {
      const t = document.createElementNS(svgNS, 'text');
      t.setAttribute('x', x); t.setAttribute('y', y);
      t.setAttribute('font-family', opts.fontFamily || 'system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial');
      t.setAttribute('font-size', (opts.fontSize || 14));
      if (opts.fontWeight) t.setAttribute('font-weight', opts.fontWeight);
      if (opts.anchor)     t.setAttribute('text-anchor', opts.anchor); // 'middle', 'start', 'end'
      if (opts.fill)       t.setAttribute('fill', opts.fill);
      t.textContent = content;
      return t;
    }

    // --- Draw header row (corner + column headers) ---
    // Corner cell
    svg.appendChild(rect(0, 0, rowHeaderW, cellH, '#ffffff', '#dddddd'));
    // Corner text
    const cornerText = text(rowHeaderW/2, textYCentered(0, cellH, headFontSize), 'Classes', {
      fontSize: headFontSize, fontWeight: 600, anchor: 'middle', fill: '#000'
    });
    svg.appendChild(cornerText);

    // Column headers
    for (let c = 1; c < cols; c++) {
      const x = rowHeaderW + (c-1)*cellW;
      const th = thead.rows[0].cells[c];
      const { primary, secondary } = parseHeaderCell(th);

      svg.appendChild(rect(x, 0, cellW, cellH, '#ffffff', '#dddddd'));

      const lineGap = Math.max(2, headFontSize * 0.2);
      const totalTextHeight = headFontSize + (secondary ? (lineGap + head2FontSize) : 0);

      // Primary centered
      const yCenter = textYCentered(0, cellH, totalTextHeight);
      const t1 = text(x + cellW/2, yCenter, primary, {
        fontSize: headFontSize, fontWeight: 600, anchor: 'middle', fill: '#000'
      });
      svg.appendChild(t1);

      if (secondary) {
        const t2 = text(x + cellW/2, yCenter + headFontSize + lineGap, `(${secondary.replace(/^\(|\)$/g,'')})`, {
          fontSize: head2FontSize, fontWeight: 500, anchor: 'middle', fill: '#6c757d'
        });
        svg.appendChild(t2);
      }
    }

    // --- Draw body (row headers + data cells) ---
    for (let r = 0; r < tbody.rows.length; r++) {
      const y = (r+1)*cellH;
      const tr = tbody.rows[r];

      // Row header cell
      const th = tr.cells[0];
      const { primary, secondary } = parseHeaderCell(th);
      svg.appendChild(rect(0, y, rowHeaderW, cellH, '#ffffff', '#dddddd'));

      const lineGap = Math.max(2, headFontSize * 0.2);
      const totalTextHeight = headFontSize + (secondary ? (lineGap + head2FontSize) : 0);
      const yCenter = textYCentered(y, cellH, totalTextHeight);

      const t1 = text(rowHeaderW/2, yCenter, primary, {
        fontSize: headFontSize, fontWeight: 600, anchor: 'middle', fill: '#000'
      });
      svg.appendChild(t1);
      if (secondary) {
        const t2 = text(rowHeaderW/2, yCenter + headFontSize + lineGap, `(${secondary.replace(/^\(|\)$/g,'')})`, {
          fontSize: head2FontSize, fontWeight: 500, anchor: 'middle', fill: '#6c757d'
        });
        svg.appendChild(t2);
      }

      // Data cells
      for (let c = 1; c < tr.cells.length; c++) {
        const x = rowHeaderW + (c-1)*cellW;
        const td = tr.cells[c];

        const bg = getComputedStyle(td).backgroundColor || 'rgb(255,255,255)';
        const fill = toRGB(bg);
        const txt = td.textContent.trim();

        svg.appendChild(rect(x, y, cellW, cellH, fill, '#dddddd'));

        if (txt) {
          const yNum = textYCentered(y, cellH, numFontSize);
          const t = text(x + cellW/2, yNum, txt, {
            fontSize: numFontSize, fontWeight: numFontWeight, anchor: 'middle', fill: '#000'
          });
          svg.appendChild(t);
        }
      }
    }

    // Serialize
    const xml = new XMLSerializer().serializeToString(svg);
    triggerDownload(xml, makeFilename(baseName, 'svg'));
  }

  global.DataTablesSVG = { export: exportTable };
})(window);
