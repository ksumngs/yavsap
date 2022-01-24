$("[data-direction]").on("click", function(e) {
    var which_function =
        $(this).data("direction") == "vertical"
        ? tree.display.spacing_x.bind(tree.display)
        : tree.display.spacing_y.bind(tree.display);
    which_function(which_function() + +$(this).data("amount")).update();
});

function sort_nodes(asc) {
    tree.resortChildren(function(a, b) {
        return (b.height - a.height || b.value - a.value) * (asc ? 1 : -1);
    });
}

$("#sort-ascending").on("click", function(e) {
    sort_nodes(true);
    tree.display.update();
});

$("#sort-descending").on("click", function(e) {
    sort_nodes(false);
    tree.display.update();
});
