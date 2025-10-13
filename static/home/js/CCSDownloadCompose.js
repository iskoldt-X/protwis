/**
 * CCSDownloadCompose.js
 * ---------------------
 * Page-specific: compose a bitmap of the table (from DataTablesImage.render)
 * with top/left padding and the axis badges, then download as PNG/JPG/PDF.
 *
 * Depends on: html2canvas, jsPDF (UMD), DataTablesImage.render
 */
(function (global) {
  if (!global.DataTablesImage || !global.DataTablesImage.render) {
    console.error('CCSDownloadCompose: DataTablesImage.render is required.');
    return;
  }
  if (typeof html2canvas === 'undefined') {
    console.error('CCSDownloadCompose: html2canvas is required.');
    return;
  }

  function makeFilename(base, ext) {
    const d = new Date(), pad = n => String(n).padStart(2, '0');
    const stamp = `${d.getFullYear()}-${pad(d.getMonth()+1)}-${pad(d.getDate())}_` +
                  `${pad(d.getHours())}-${pad(d.getMinutes())}-${pad(d.getSeconds())}`;
    return `${base || 'ccs-heatmap'}-${stamp}.${ext}`;
  }
  function triggerDownload(dataUrl, filename) {
    const a = document.createElement('a');
    a.href = dataUrl; a.download = filename; document.body.appendChild(a); a.click(); a.remove();
  }
  async function captureBadgeCanvas(el, scale, { fixVertical=false } = {}) {
    if (!el) return null;

    // Clone off-screen to avoid mutating the live DOM
    const clone = el.cloneNode(true);
    clone.style.position = 'fixed';
    clone.style.left = '-10000px';
    clone.style.top = '0';
    clone.style.margin = '0';

    // When capturing the Y badge, avoid writing-mode issues in html2canvas:
    // render it as horizontally-laid-out text rotated -90deg instead.
    if (fixVertical) {
      // Neutralize original vertical styling
      clone.style.writingMode = 'horizontal-tb';
      clone.style.transform = 'rotate(-90deg)';       // rotate the whole pill
      clone.style.transformOrigin = 'left top';
      clone.style.display = 'inline-block';
    }

    document.body.appendChild(clone);

    // Let layout settle before snapshot
    await new Promise(r => requestAnimationFrame(r));

    const canvas = await html2canvas(clone, {
      backgroundColor: null,
      useCORS: true,
      scale: scale || Math.max(2, (window.devicePixelRatio || 1) * 2),
      scrollX: 0,
      scrollY: 0
    });

    document.body.removeChild(clone);
    return canvas;
  }

  async function composeAndDownload(format, userOpts) {
    const opts = Object.assign({
      tableSelector: '#ccsTable',
      axisWrapSelector: '.ccs-axis-wrap',
      padTop: 60,
      padLeft: 60,
      gapX: 8,     // gap between table top and X-badge (visual)
      gapY: 10,    // gap between table left and Y-badge (visual)
      scale: Math.max(2, (window.devicePixelRatio || 1) * 2),
      quality: 0.95,
      baseName: 'ccs-heatmap',

      // backgrounds:
      tableBackground: '#ffffff',   // the table itself stays opaque with this color
      background: null,             // outer padding: null/'transparent' => transparent; or a CSS color
      pdfBackground: null           // optional PDF page fill color
    }, userOpts || {});

    const ext = String(format || 'png').toLowerCase();
    const wantsTransparentOuter =
      opts.background === null || String(opts.background).toLowerCase() === 'transparent';

    // 1) Render the table as a canvas (opaque background so table is not transparent)
    const tableCanvas = await window.DataTablesImage.render(
      opts.tableSelector,
      { scale: opts.scale, background: opts.tableBackground }
    );

    // 2) Capture the badges as canvases (they naturally have transparent bg in our capture)
    const axisWrap = document.querySelector(opts.axisWrapSelector);
    const xBadgeEl = axisWrap && axisWrap.querySelector('.axis-label-x .ccs-axis-badge');
    const yBadgeEl = axisWrap && axisWrap.querySelector('.axis-label-y .ccs-axis-badge');

    // NOTE: pass { fixVertical:true } for the Y badge to avoid writing-mode issues in html2canvas
    const [xBadgeCanvas, yBadgeCanvas] = await Promise.all([
      captureBadgeCanvas(xBadgeEl, opts.scale),
      captureBadgeCanvas(yBadgeEl, opts.scale, { fixVertical: true })
    ]);

    // 3) Compute padding so badges fit if they are larger than the requested pad
    const padTop  = Math.max(opts.padTop,  (xBadgeCanvas ? xBadgeCanvas.height + opts.gapX : 0));
    const padLeft = Math.max(opts.padLeft, (yBadgeCanvas ? yBadgeCanvas.width  + opts.gapY : 0));

    // 4) Compose onto a new canvas
    const W = padLeft + tableCanvas.width;
    const H = padTop  + tableCanvas.height;

    const out = document.createElement('canvas');
    out.width = W; out.height = H;
    const ctx = out.getContext('2d');

    // Fill outer background only if requested; otherwise keep outer fully transparent
    if (!wantsTransparentOuter) {
      ctx.fillStyle = opts.background || '#ffffff';
      ctx.fillRect(0, 0, W, H);
    } else {
      ctx.clearRect(0, 0, W, H);
    }

    // Draw the X badge centered above the table
    if (xBadgeCanvas) {
      const bx = padLeft + Math.round((tableCanvas.width - xBadgeCanvas.width) / 2);
      const by = Math.max(0, padTop - opts.gapX - xBadgeCanvas.height);
      ctx.drawImage(xBadgeCanvas, bx, by);
    }

    // Draw the Y badge to the left, vertically centered to the table
    if (yBadgeCanvas) {
      const bx = Math.max(0, padLeft - opts.gapY - yBadgeCanvas.width);
      const by = padTop + Math.round((tableCanvas.height - yBadgeCanvas.height) / 2);
      ctx.drawImage(yBadgeCanvas, bx, by);
    }

    // Draw the table (already opaque)
    ctx.drawImage(tableCanvas, padLeft, padTop);

    // 5) Save
    const filename = (function() {
      const d = new Date(), pad = n => String(n).padStart(2, '0');
      const stamp = `${d.getFullYear()}-${pad(d.getMonth()+1)}-${pad(d.getDate())}_` +
                    `${pad(d.getHours())}-${pad(d.getMinutes())}-${pad(d.getSeconds())}`;
      return `${opts.baseName}-${stamp}.${ext}`;
    })();

    if (ext === 'png') {
      // PNG preserves transparency in the outer padding if wantsTransparentOuter === true
      const url = out.toDataURL('image/png');
      const a = document.createElement('a'); a.href = url; a.download = filename; a.click();
    } else if (ext === 'jpg' || ext === 'jpeg') {
      // JPEG has no alpha; if outer is transparent, flatten onto white (or a chosen color)
      if (wantsTransparentOuter) {
        const flat = document.createElement('canvas');
        flat.width = W; flat.height = H;
        const fctx = flat.getContext('2d');
        fctx.fillStyle = '#ffffff';
        fctx.fillRect(0, 0, W, H);
        fctx.drawImage(out, 0, 0);
        const url = flat.toDataURL('image/jpeg', opts.quality);
        const a = document.createElement('a'); a.href = url; a.download = filename; a.click();
      } else {
        const url = out.toDataURL('image/jpeg', opts.quality);
        const a = document.createElement('a'); a.href = url; a.download = filename; a.click();
      }
    } else if (ext === 'pdf') {
      if (typeof window.jspdf === 'undefined') {
        console.error('CCSDownloadCompose: jsPDF is required for PDF output.');
        return;
      }
      const { jsPDF } = window.jspdf;
      const pdf = new jsPDF({ orientation: (W >= H) ? 'l' : 'p', unit: 'pt', format: [W, H] });

      // Optional PDF page background fill (useful if you want colored page behind transparent PNG)
      if (opts.pdfBackground) {
        pdf.setFillColor(opts.pdfBackground);
        pdf.rect(0, 0, W, H, 'F');
      }

      pdf.addImage(out.toDataURL('image/png'), 'PNG', 0, 0, W, H);
      pdf.save(filename);
    } else {
      console.error('CCSDownloadCompose: unsupported format:', ext);
    }
  }


  // public API
  global.CCSDownload = {
    withBadges: composeAndDownload
  };
})(window);
