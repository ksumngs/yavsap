$("[data-direction]").on("click", function(e) {
    var which_function =
      $(this).data("direction") == "vertical"
        ? tree.display.spacing_x.bind(tree.display)
        : tree.display.spacing_y.bind(tree.display);
    which_function(which_function() + +$(this).data("amount")).update();
  });
