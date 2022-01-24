function zoom(amount) {
    xzoom = tree.display.spacing_x.bind(tree.display);
    yzoom = tree.display.spacing_y.bind(tree.display);
    xzoom(xzoom() + amount).update();
    yzoom(yzoom() + amount).update();
}

function zoomin(e) {
    zoom(1);
}

function zoomout(e) {
    zoom(-1);
}
