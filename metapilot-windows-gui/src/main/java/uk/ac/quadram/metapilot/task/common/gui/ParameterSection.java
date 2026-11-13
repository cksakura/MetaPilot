package uk.ac.quadram.metapilot.task.common.gui;

import javafx.geometry.Insets;
import javafx.scene.Node;
import javafx.scene.control.ScrollPane;
import javafx.scene.control.TitledPane;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.Priority;
import javafx.scene.layout.StackPane;
import javafx.scene.layout.VBox;

public final class ParameterSection {

    private ParameterSection() {
    }

    public static ScrollPane createViewport(GridPane grid) {
        ScrollPane viewport = new ScrollPane(grid);
        viewport.setFitToWidth(true);
        viewport.setFitToHeight(true);
        viewport.setStyle("-fx-background-insets: 0; -fx-padding: 0;");
        VBox.setVgrow(viewport, Priority.ALWAYS);
        return viewport;
    }

    public static TitledPane create(String title, Node content) {
        StackPane body = new StackPane(content);
        body.setPadding(new Insets(10));
        TitledPane section = new TitledPane(title, body);
        section.getStyleClass().add("parameter-section");
        section.setCollapsible(false);
        section.setMaxHeight(Double.MAX_VALUE);
        return section;
    }
}
