/* Pause the drifting sequences behind the landing hero while the hero is
 * off screen, so the page does no animation work nobody can see.
 *
 * The rows themselves are plain HTML and the motion is pure CSS (see
 * landing.css), so the hero looks complete without this file. Under
 * prefers-reduced-motion the CSS turns the animation off altogether.
 */
(function () {
  "use strict";

  function watch(hero) {
    if (!("IntersectionObserver" in window)) {
      return;
    }
    var observer = new IntersectionObserver(function (entries) {
      for (var i = 0; i < entries.length; i++) {
        hero.classList.toggle("is-offscreen", !entries[i].isIntersecting);
      }
    });
    observer.observe(hero);
  }

  document.addEventListener("DOMContentLoaded", function () {
    var hero = document.querySelector(".rot-hero");
    if (hero) {
      watch(hero);
    }
  });
})();
