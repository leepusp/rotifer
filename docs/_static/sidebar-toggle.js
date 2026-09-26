/* Collapse and restore both sidebars.
 *
 * Loaded from <head> (see _templates/partials/extra-head.html) so a
 * remembered state is on <html> before the first paint. Everything visual
 * lives in _static/theme.css, keyed off the classes this file sets.
 *
 * The listeners are bound to the buttons themselves rather than delegated
 * from the document, and that is not a style choice. Shibuya wires its
 * drawer buttons with a helper that registers a click handler ON the aside
 * which calls stopPropagation(), so no click inside a sidebar ever reaches
 * the document. A delegated handler silently never fires.
 *
 * Each panel has a width below which the theme turns it into an overlay
 * drawer with its own close button. Below that width the stylesheet makes
 * the control inert and this file ignores it, so the two mechanisms never
 * fight over the same panel.
 */
(function () {
  "use strict";

  var PANELS = [
    {
      /* Left: navigation. In flow from 768px up. */
      hook: "js-lside-toggle",
      klass: "lside-collapsed",
      key: "rotifer-lside-collapsed",
      minWidth: 768
    },
    {
      /* Right: "On this page". In flow from 1280px up. */
      hook: "js-rside-toggle",
      klass: "rside-collapsed",
      key: "rotifer-localtoc-collapsed",
      minWidth: 1280
    }
  ];

  /* Long enough to read as one movement, short enough not to be in the
   * way. Must match --rot-transition-panel in theme.css. */
  var ANIMATING = "sidebars-animating";
  var DURATION = 500;

  var root = document.documentElement;
  var timer = null;

  /* Private windows and blocked site data make localStorage throw, and the
   * sidebars have to work anyway, so every access is guarded and expanded
   * is the fallback. */
  function remembered(panel) {
    try {
      return window.localStorage.getItem(panel.key) === "1";
    } catch (err) {
      return false;
    }
  }

  function remember(panel, collapsed) {
    try {
      window.localStorage.setItem(panel.key, collapsed ? "1" : "0");
    } catch (err) {
      /* Nothing to do: the panel still works for this page view. */
    }
  }

  function buttons(panel) {
    return document.querySelectorAll("." + panel.hook);
  }

  function apply(panel, collapsed) {
    root.classList.toggle(panel.klass, collapsed);
    var found = buttons(panel);
    for (var i = 0; i < found.length; i++) {
      found[i].setAttribute("aria-expanded", collapsed ? "false" : "true");
    }
  }

  /* The width transition is switched on only for the length of a toggle.
   * Left on permanently it would also animate every window resize, which
   * reads as lag rather than as motion. */
  function animate() {
    root.classList.add(ANIMATING);
    if (timer !== null) {
      window.clearTimeout(timer);
    }
    timer = window.setTimeout(function () {
      root.classList.remove(ANIMATING);
      timer = null;
    }, DURATION);
  }

  function toggle(panel) {
    if (window.innerWidth < panel.minWidth) {
      return;
    }
    var collapsed = !root.classList.contains(panel.klass);
    animate();
    apply(panel, collapsed);
    remember(panel, collapsed);
  }

  function bind(panel) {
    var found = buttons(panel);
    for (var i = 0; i < found.length; i++) {
      found[i].addEventListener("click", function () {
        toggle(panel);
      });
    }
    /* The buttons exist now, so give them their aria state. */
    apply(panel, root.classList.contains(panel.klass));
  }

  for (var i = 0; i < PANELS.length; i++) {
    apply(PANELS[i], remembered(PANELS[i]));
  }

  document.addEventListener("DOMContentLoaded", function () {
    for (var j = 0; j < PANELS.length; j++) {
      bind(PANELS[j]);
    }
  });
})();
