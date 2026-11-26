/* --------------------------------------------------------------------------
   NADIA — pkgdown site
   --------------------------------------------------------------------------
   Figure lightbox for the articles. A vignette marks a figure as enlargeable by
   wrapping it in a `.lightbox-figure` div; the class is inert in the BiocStyle
   build of the vignettes and only does something here. Each vignette used to
   carry its own copy of this script (Bioconductor review, #186).

   Clicking the image, or pressing Enter or Space on it, opens it at the size of
   the window; Esc, the close button or a click outside the image closes it.
   The hint under each figure is added here too, so that it never appears
   where there is no lightbox.
   -------------------------------------------------------------------------- */

(function () {
  function initialiseLightbox() {
    var images = document.querySelectorAll(".lightbox-figure img");
    if (!images.length) return;

    var lightbox = document.createElement("div");
    lightbox.id = "nadia-lightbox";
    lightbox.setAttribute("role", "dialog");
    lightbox.setAttribute("aria-modal", "true");
    lightbox.setAttribute("aria-label", "Enlarged figure");
    lightbox.setAttribute("aria-hidden", "true");

    var enlarged = document.createElement("img");
    enlarged.alt = "";
    var closeButton = document.createElement("button");
    closeButton.type = "button";
    closeButton.setAttribute("aria-label", "Close enlarged figure");
    closeButton.textContent = "×";

    lightbox.appendChild(enlarged);
    lightbox.appendChild(closeButton);
    document.body.appendChild(lightbox);

    var opener = null;

    function openLightbox(source) {
      opener = source;
      enlarged.src = source.currentSrc || source.src;
      enlarged.alt = source.alt || "Enlarged figure";
      lightbox.classList.add("is-open");
      lightbox.setAttribute("aria-hidden", "false");
      document.body.classList.add("nadia-lightbox-open");
      closeButton.focus();
    }

    function closeLightbox() {
      if (!lightbox.classList.contains("is-open")) return;
      lightbox.classList.remove("is-open");
      lightbox.setAttribute("aria-hidden", "true");
      document.body.classList.remove("nadia-lightbox-open");
      if (opener) opener.focus();
    }

    images.forEach(function (source) {
      source.setAttribute("tabindex", "0");
      source.setAttribute("title", "Click to enlarge");
      source.addEventListener("click", function () {
        openLightbox(source);
      });
      source.addEventListener("keydown", function (event) {
        if (event.key === "Enter" || event.key === " ") {
          event.preventDefault();
          openLightbox(source);
        }
      });
    });

    document.querySelectorAll(".lightbox-figure").forEach(function (figure) {
      // A wrapper class ending in "-table" holds an image of a table.
      var what = /(^|\s)[\w-]+-table(\s|$)/.test(figure.className) ? "table" : "figure";
      var hint = document.createElement("p");
      hint.className = "figure-zoom-hint";
      hint.textContent = "Click the " + what + " to enlarge it; press Esc or " +
        "click outside the image to close it.";
      figure.parentNode.insertBefore(hint, figure.nextSibling);
    });

    enlarged.addEventListener("click", function (event) {
      event.stopPropagation();
    });
    closeButton.addEventListener("click", closeLightbox);
    lightbox.addEventListener("click", closeLightbox);
    document.addEventListener("keydown", function (event) {
      if (event.key === "Escape") closeLightbox();
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", initialiseLightbox);
  } else {
    initialiseLightbox();
  }
})();
