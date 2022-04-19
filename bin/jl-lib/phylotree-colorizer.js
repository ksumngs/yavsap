tree = new phylotree.phylotree(newick);
selection_set = !tree.get_parsed_tags().length
  ? tree.get_parsed_tags()
  : ["Foreground"];
nodeColorizer = function (element, data) {
  try {
    var count_class = 0;
    selection_set.forEach(function (d, i) {
      if (data[d]) {
        count_class++;
        element.style("fill", color_scale(i), "important");
      }
    });
    if (count_class == 0) {
      element.style("fill", null);
    }
  } catch (e) {}
};
edgeColorizer = function (element, data) {
  try {
    var count_class = 0;

    selection_set.forEach(function (d, i) {
      if (data[d]) {
        count_class++;
        element.style("stroke", color_scale(i), "important");
      }
    });

    if (count_class == 0) {
      element.style("stroke", null).classes("branch-multiple", false);
    } else {
      element.classed("branch-multiple", true);
    }
  } catch (e) {}
};
colorNodesByName = function (element, data) {
  nodeColorizer(element, data);
  var m = data.data.name.split("_");
  element.style("stroke", color_scale(m[0]));
};
colorEdgesByTarget = function (element, data) {
  edgeColorizer(element, data);
  var m = data.target.data.name.split("_");
  element.style("stroke", color_scale(m[0]));
};
rendered_tree = tree.render({
  container: "#phylotree",
  "node-styler": colorNodesByName,
  //'edge-styler': colorEdgesByTarget
});
$(tree.display.container).empty();
$(tree.display.container).html(tree.display.show());
