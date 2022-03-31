/*
  phylotree-utils.js

  Extra functions that make Phylotree.js more useful.
  Ported from Phylotree.js's example document.
  Copyright (c) 2016 iGEM/UCSD evolutionary biology and bioinformatics group
  https://github.com/veg/phylotree.js/blob/93fdebb81503f83b3fffe0a56ad3c02c64535fea/index.html

  Used under MIT License

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
*/

$("[data-direction]").on("click", function (e) {
  var which_function =
    $(this).data("direction") == "vertical"
      ? tree.display.spacing_x.bind(tree.display)
      : tree.display.spacing_y.bind(tree.display);
  which_function(which_function() + +$(this).data("amount")).update();
});

function sort_nodes(asc) {
  tree.resortChildren(function (a, b) {
    return (b.height - a.height || b.value - a.value) * (asc ? 1 : -1);
  });
}

$("#sort-ascending").on("click", function (e) {
  sort_nodes(true);
  tree.display.update();
});

$("#sort-descending").on("click", function (e) {
  sort_nodes(false);
  tree.display.update();
});

var datamonkey_save_image = function (type, container) {
  var prefix = {
    xmlns: "http://www.w3.org/2000/xmlns/",
    xlink: "http://www.w3.org/1999/xlink",
    svg: "http://www.w3.org/2000/svg",
  };

  var svg = $(container).find("svg")[0];
  if (!svg) {
    svg = $(container)[0];
  }

  svg.setAttribute("version", "1.1");

  var defsEl = document.createElement("defs");
  svg.insertBefore(defsEl, svg.firstChild);

  var styleEl = document.createElement("style");
  defsEl.appendChild(styleEl);
  styleEl.setAttribute("type", "text/css");

  // removing attributes so they aren't doubled up
  svg.removeAttribute("xmlns");
  svg.removeAttribute("xlink");

  // These are needed for the svg
  if (!svg.hasAttributeNS(prefix.xmlns, "xmlns")) {
    svg.setAttributeNS(prefix.xmlns, "xmlns", prefix.svg);
  }

  if (!svg.hasAttributeNS(prefix.xmlns, "xmlns:xlink")) {
    svg.setAttributeNS(prefix.xmlns, "xmlns:xlink", prefix.xlink);
  }

  // This is a deviation from the original Phylotree-utils
  // We know that the 2nd stylesheet contains all the svg styles,
  // so import them
  styles = "";
  styleRules = [...document.styleSheets[1].cssRules];
  styleRules.forEach(function (style) {
    styles = styles + style.cssText + "\n";
  });
  svg.getElementsByTagName("style")[0].innerText = styles;

  var source = new XMLSerializer().serializeToString(svg);
  var doctype = new XMLSerializer().serializeToString(
    document.implementation.createDocumentType(
      "svg",
      "-//W3C//DTD SVG 1.1//EN",
      "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd"
    )
  );
  var to_download = [doctype + source];
  var image_string =
    "data:image/svg+xml;base66," + encodeURIComponent(to_download);

  if (navigator.msSaveBlob) {
    // IE10
    download(image_string, "image.svg", "image/svg+xml");
  } else if (type == "png") {
    b64toBlob(
      image_string,
      function (blob) {
        var url = window.URL.createObjectURL(blob);
        var pom = document.createElement("a");
        pom.setAttribute("download", "image.png");
        pom.setAttribute("href", url);
        $("body").append(pom);
        pom.click();
        pom.remove();
      },
      function (error) {
        console.log(error); // eslint-disable-line
      }
    );
  } else {
    var pom = document.createElement("a");
    pom.setAttribute("download", "image.svg");
    pom.setAttribute("href", image_string);
    $("body").append(pom);
    pom.click();
    pom.remove();
  }
};

$("#save-image").on("click", function (e) {
  datamonkey_save_image("svg", "#phylotree");
});

$(".phylotree-layout-mode").on("click", function (e) {
  if (tree.display.radial() != ($(this).data("mode") == "radial")) {
    $(".phylotree-layout-mode").toggleClass("active");
    tree.display.radial(!tree.display.radial()).update();
  }
});

$(".phylotree-align-toggler").on("click", function (e) {
  var button_align = $(this).data("align");
  var tree_align = tree.display.options.alignTips;

  if (tree_align != button_align) {
    tree.display.alignTips(button_align == "right");
    $(".phylotree-align-toggler").toggleClass("active");
    tree.display.update();
  }
});
color_scale = d3.scaleOrdinal(d3.schemeCategory10);
