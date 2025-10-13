/**
 * DataTablesImage.js
 * ------------------
 * Export or render a single DOM element (e.g., a <table>) with a tight crop.
 *
 * Usage:
 *   // Download directly
 *   DataTablesImage.export('#ccsTable','ccs-table',{ format:'png', scale:3 })
 *
 *   // Get a canvas without downloading (for custom composition)
 *   const canvas = await DataTablesImage.render('#ccsTable', { scale:3, background:'#fff' });
 *
 * Depends on: html2canvas, jsPDF (UMD for PDF only)
 */

(function (global) {
  if (typeof html2canvas === 'undefined') {
    console.error('DataTablesImage.js: html2canvas is required.');
    return;
  }

  function triggerDownload(dataUrl, filename) {
    const a = document.createElement('a');
    a.href = dataUrl;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    a.remove();
  }

  function makeFilename(base, ext) {
    const d = new Date(), pad = n => String(n).padStart(2, '0');
    const stamp = `${d.getFullYear()}-${pad(d.getMonth()+1)}-${pad(d.getDate())}_` +
                  `${pad(d.getHours())}-${pad(d.getMinutes())}-${pad(d.getSeconds())}`;
    return `${base || 'table'}-${stamp}.${ext}`;
  }

  // Temporarily relax the immediate parent's overflow so html2canvas can paint full content
  function loosenOverflow(el) {
    const parent = el && el.parentElement;
    if (!parent) return () => {};
    const prev = {
      overflow: parent.style.overflow,
      overflowX: parent.style.overflowX,
      overflowY: parent.style.overflowY,
      maxWidth: parent.style.maxWidth,
      maxHeight: parent.style.maxHeight
    };
    parent.style.overflow = 'visible';
    parent.style.overflowX = 'visible';
    parent.style.overflowY = 'visible';
    parent.style.maxWidth = 'none';
    parent.style.maxHeight = 'none';
    return () => {
      parent.style.overflow   = prev.overflow;
      parent.style.overflowX  = prev.overflowX;
      parent.style.overflowY  = prev.overflowY;
      parent.style.maxWidth   = prev.maxWidth;
      parent.style.maxHeight  = prev.maxHeight;
    };
  }

  /**
   * Render only (no download). Returns a canvas tightly cropped to the element.
   * @param {string|HTMLElement} selector
   * @param {{scale?:number, background?:string}} opts
   * @returns {Promise<HTMLCanvasElement>}
   */
  async function renderElement(selector, opts) {
    const options = Object.assign(
      { scale: Math.max(2, (window.devicePixelRatio || 1) * 2), background: '#ffffff' },
      opts || {}
    );

    const el = (typeof selector === 'string') ? document.querySelector(selector) : selector;
    if (!el) throw new Error('DataTablesImage.render: selector not found: ' + selector);

    const restore = loosenOverflow(el);
    await new Promise(r => requestAnimationFrame(r)); // let layout settle

    const canvas = await html2canvas(el, {
      backgroundColor: options.background,
      useCORS: true,
      scale: options.scale,
      scrollX: 0,
      scrollY: 0
    });

    restore();
    return canvas;
  }

  /**
   * Render and download as PNG/JPG/PDF.
   * @param {string|HTMLElement} selector
   * @param {string} baseName
   * @param {{format?:'png'|'jpg'|'jpeg'|'pdf', quality?:number, scale?:number}} opts
   */
  async function exportElement(selector, baseName, opts) {
    const options = Object.assign(
      { format: 'png', quality: 0.92, scale: Math.max(2, (window.devicePixelRatio || 1) * 2) },
      opts || {}
    );
    const fmt = String(options.format || 'png').toLowerCase();

    const canvas = await renderElement(selector, { scale: options.scale, background: '#ffffff' });
    const filename = makeFilename(baseName, fmt);

    if (fmt === 'png') {
      triggerDownload(canvas.toDataURL('image/png'), filename);
    } else if (fmt === 'jpg' || fmt === 'jpeg') {
      triggerDownload(canvas.toDataURL('image/jpeg', options.quality), filename);
    } else if (fmt === 'pdf') {
      if (typeof window.jspdf === 'undefined') {
        console.error('DataTablesImage.export: jsPDF is required for PDF output.');
        return;
      }
      const { jsPDF } = window.jspdf;
      const w = canvas.width, h = canvas.height;
      const pdf = new jsPDF({ orientation: (w >= h) ? 'l' : 'p', unit: 'pt', format: [w, h] });
      pdf.addImage(canvas.toDataURL('image/png'), 'PNG', 0, 0, w, h);
      pdf.save(filename);
    } else {
      console.error('DataTablesImage.export: unsupported format:', fmt);
    }
  }

  // Public API
  global.DataTablesImage = {
    render: renderElement,
    export: exportElement
  };
})(window);
