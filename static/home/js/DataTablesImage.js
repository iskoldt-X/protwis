/**
 * DataTablesImage.js
 * ------------------
 * Export or render a single DOM element (e.g., a <table>) with a tight crop.
 *
 * Usage:
 *   DataTablesImage.export('#ccsTable','ccs-table',{ format:'png', scale:3, filename:'MyFile.png' })
 *   const canvas = await DataTablesImage.render('#ccsTable', { scale:3 });
 *
 * Depends on: html2canvas
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

  function basenameFromSelector(selector, fallback) {
    if (typeof selector === 'string') {
      const m = selector.match(/#([\w-]+)/);
      if (m && m[1]) return m[1];
    }
    return fallback || 'table';
  }

  function stamped(base, ext) {
    const d = new Date(), pad = n => String(n).padStart(2, '0');
    const stamp = `${d.getFullYear()}-${pad(d.getMonth()+1)}-${pad(d.getDate())}_` +
                  `${pad(d.getHours())}-${pad(d.getMinutes())}-${pad(d.getSeconds())}`;
    return `${base}-${stamp}.${ext}`;
  }

  function resolveFilename(selector, baseName, ext, explicit) {
    if (explicit) return explicit; // honor explicit filename
    const base = (baseName && String(baseName).trim()) || basenameFromSelector(selector, 'table');
    return stamped(base, ext);
  }

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

  async function renderElement(selector, opts) {
    const options = Object.assign(
      { scale: Math.max(2, (window.devicePixelRatio || 1) * 2), background: '#ffffff' },
      opts || {}
    );
    const el = (typeof selector === 'string') ? document.querySelector(selector) : selector;
    if (!el) throw new Error('DataTablesImage.render: selector not found: ' + selector);

    const restore = loosenOverflow(el);
    await new Promise(r => requestAnimationFrame(r));
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

  async function exportElement(selector, baseName, opts) {
    const options = Object.assign(
      { format: 'png', quality: 0.92, scale: Math.max(2, (window.devicePixelRatio || 1) * 2), background:'#ffffff', filename:null },
      opts || {}
    );
    const fmt = String(options.format || 'png').toLowerCase();
    const canvas = await renderElement(selector, { scale: options.scale, background: options.background });
    const filename = resolveFilename(selector, baseName, fmt, options.filename);

    if (fmt === 'png') {
      triggerDownload(canvas.toDataURL('image/png'), filename);
    } else if (fmt === 'jpg' || fmt === 'jpeg') {
      triggerDownload(canvas.toDataURL('image/jpeg', options.quality), filename);
    } else {
      console.error('DataTablesImage.export: unsupported format (use png/jpg):', fmt);
    }
  }

  global.DataTablesImage = {
    render: renderElement,
    export: exportElement
  };
})(window);
