help_modal = html_div(
    html_div(
        html_div(
            header(
                h5("Help"; class="modal-title"),
                button(
                    "";
                    type="button",
                    class="btn-close",
                    data_bs_dismiss="modal",
                    aria_label="Close",
                );
                class="modal-header",
            );
            class="modal-content",
        ),
        ;
        class="modal-dialog modal-dialog-centered",
    ),
    html_div("Lorum ipsum text..."; class="modal-body");
    id="help-dialog",
    class="modal fade",
    tabindex="-1",
)

phylogenetic_section = section(
    h2("phylotree.js"),
    html_div(
        nav(
            html_div(
                button(
                    i(""; class="fas fa-arrows-alt-v");
                    class="btn btn-primary",
                    type="button",
                    data_direction="vertical",
                    data_amount="1",
                    data_bs_toggle="tooltip",
                    title="Expand vertical spacing",
                ),
                button(
                    i(
                        "";
                        class="fas fa-compress-alt fa-rotate-by",
                        style="--fa-rotate-angle: 135deg;",
                    );
                    class="btn btn-primary",
                    type="button",
                    data_direction="vertical",
                    data_amount="-1",
                    data_bs_toggle="tooltip",
                    title="Compress vertical spacing",
                ),
                button(
                    i(""; class="fas fa-arrows-alt-h");
                    class="btn btn-primary",
                    type="button",
                    data_direction="horizontal",
                    data_amount="1",
                    data_bs_toggle="tooltip",
                    title="Expand horizontal spacing",
                ),
                button(
                    i(
                        "";
                        class="fas fa-compress-alt fa-rotate-by",
                        style="--fa-rotate-angle: 45deg;",
                    );
                    class="btn btn-primary",
                    type="button",
                    data_direction="horizontal",
                    data_amount="-1",
                    data_bs_toggle="tooltip",
                    title="Compress horizontal spacing",
                );
                class="btn-group mx-1",
                role="group",
            ),
            html_div(
                button(
                    i(""; class="fas fa-sort-amount-down-alt");
                    id="sort-ascending",
                    class="btn btn-primary",
                    type="button",
                    data_bs_toggle="tooltip",
                    title="Sort deepest clades to the bottom",
                ),
                button(
                    i(""; class="fas fa-sort-amount-up-alt");
                    id="sort-descending",
                    class="btn btn-primary",
                    type="button",
                    data_bs_toggle="tooltip",
                    title="Sort deepest clades to the top",
                );
                class="btn-group mx-1",
                role="group",
            ),
            html_div(
                button(
                    i(""; class="fas fa-file-image");
                    id="save-image",
                    class="btn btn-primary",
                    type="button",
                    data_bs_toggle="tooltip",
                    title="Save image",
                );
                class="btn-group mx-1",
            ),
            html_div(
                button(
                    i(""; class="fas fa-align-left");
                    class="btn btn-primary phylotree-layout-mode active",
                    type="button",
                    title="Linear",
                    data_mode="linear",
                ),
                button(
                    i(""; class="fas fa-circle-notch");
                    class="btn btn-primary phylotree-layout-mode",
                    type="button",
                    title="Radial",
                    data_mode="radial",
                );
                class="btn-group mx-1",
            ),
            html_div(
                button(
                    i(""; class="fas fa-outdent");
                    class="btn btn-primary phylotree-align-toggler active",
                    type="button",
                    title="Align left",
                    data_align="left",
                ),
                html_div(
                    button(
                        i(""; class="fas fa-indent");
                        class="btn btn-primary phylotree-align-toggler",
                        type="button",
                        title="Align right",
                        data_align="right",
                    ),
                );
                class="btn-group mx-1",
                role="group",
            );
            class="col-12",
        ),
        html_div(html_div(""; id="phylotree", class="my-2 p-3"); class="col-12");
        class="row",
    );
    id="phylogenetics",
    class="container-fluid min-vh-100",
)
