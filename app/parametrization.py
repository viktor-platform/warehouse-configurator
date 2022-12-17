from viktor.geometry import GeoPoint, GeoPolygon
from viktor.parametrization import (
    GeoPointField,
    GeoPolygonField,
    IntegerField,
    LineBreak,
    Lookup,
    NumberField,
    OptionField,
    OptionListElement,
    Parametrization,
    Section,
    Step,
    Tab,
    Text,
    ToggleButton,
)

from .database import profiles_options

DEFAULT_CORNER_LOCATION = GeoPoint(52.27708034013321, 4.749773456312529)

DEFAULT_TERRAIN = GeoPolygon(
    GeoPoint(52.277114025504986, 4.748731333270152),
    GeoPoint(52.27664483409234, 4.7496200872609435),
    GeoPoint(52.27695522281373, 4.7500998571143676),
    GeoPoint(52.27740516229457, 4.749163912646231),
)

PURLIN_SPACING_OPTIONS = [
    OptionListElement(label="1/2", value=2.0),
    OptionListElement(label="1", value=1.0),
    OptionListElement(label="2", value=0.5),
]

beam_type = [
    "HEA",
    "HEB",
    "IPE",
    "SHS",
]


class WarehouseParametrization(Parametrization):

    location_and_building = Step(
        "Location and Buiding", views=["get_map_view", "visualize_exterior"]
    )
    location_and_building.text_01 = Text(
        """
# Welcome to the Warehouse Configurator app!

This app demonstrates how a user could design a warehouse with an adjacent office in a simplified an parametrized 
manner.

## Location

Start of by selecting the terrain where the building should be situated. Select then the location of the warehouse, 
as well as its orientation.
    """
    )
    location_and_building.poly = GeoPolygonField(
        "Terrain Polygon", default=DEFAULT_TERRAIN
    )
    location_and_building.start = GeoPointField(
        "Building corner", default=DEFAULT_CORNER_LOCATION
    )
    location_and_building.rotate = NumberField("Rotate", suffix="°", default=230)

    location_and_building.building_text = Text(
        """
## Building

For this design, the warehouse and office share a "wall" with the same length.

Click on the "Building" tab on the right to see a simplified rendered model of the warehouse design.
    """
    )
    location_and_building.building_y_dimension = NumberField(
        "Office/warehouse length", min=20, default=20, step=10, suffix="m"
    )

    location_and_building.office_text = Text("### Office")
    location_and_building.office_x_dimension = NumberField(
        "Office width", min=10, default=20, step=10, suffix="m"
    )
    location_and_building.num_office_floors = IntegerField(
        "Number of floors", min=3, default=3, step=1
    )

    location_and_building.warehouse_text = Text("### Warehouse")
    location_and_building.warehouse_x_dimension = NumberField(
        "Warehouse width", min=20, default=20, step=10, suffix="m"
    )
    location_and_building.num_warehouse_floors = IntegerField(
        "Free height",
        min=2,
        default=2,
        suffix="floors",
        description="The vertical free space for the warehouse is defined in the equivalent number of floors",
    )

    structure = Step("Steel Frame", views="visualize_structure")
    structure.general = Tab("General")

    structure.general.office_col = Section("Office columns")

    structure.general.office_col.dist_length = NumberField(
        "Max column spacing, length", min=1, default=7, step=0.5, suffix="m"
    )
    structure.general.office_col.dist_width = NumberField(
        "M column spacing, width", min=1, default=7, step=0.5, suffix="m"
    )
    structure.general.office_col.lb2 = LineBreak()
    structure.general.office_col.col_profile = OptionField(
        "Column profile", options=profiles_options, default="SHS 300x300 x 10"
    )

    structure.general.truss = Section("Warehouse trusses")
    structure.general.truss.text = Text(
        "The trusses are used to support the roof for long spans."
    )
    structure.general.truss.max_truss_spacing = NumberField(
        "Max. truss spacing",
        min=5,
        default=15,
        step=5,
        suffix="m",
        description="The maximum distance between two trusses.",
    )
    structure.general.truss.lb1 = LineBreak()
    structure.general.truss.truss_height = NumberField(
        "Truss height", min=1, default=1.0, step=0.2, suffix="m"
    )
    structure.general.truss.lb2 = LineBreak()
    structure.general.truss.profile_chord = OptionField(
        "Truss chord profile",
        options=profiles_options,
        default="SHS 100x100 x 4",
        description="Profile of the horizontal beam.",
    )
    structure.general.truss.profile_web = OptionField(
        "Truss web profile",
        options=profiles_options,
        default="SHS 50x50 x 4",
        description="Profile of the diagonal beam.",
    )
    structure.general.truss.profile_vertical = OptionField(
        "Truss vertical profile", options=profiles_options, default="SHS 50x50 x 4"
    )
    structure.general.truss.lb3 = LineBreak()
    structure.general.truss.custom_panels = ToggleButton("Custom panels")
    structure.general.truss.truss_panels = IntegerField(
        "Panels per truss",
        min=2,
        step=2,
        default=2,
        visible=Lookup("structure.general.truss.sqr_panels"),
    )

    structure.general.column = Section("Warehouse columns")
    structure.general.column.text = Text(
        "Whenever the trusses cannot hold the weight of the roof, "
        "extra columns can be added to the structure."
    )
    structure.general.column.num_columns = IntegerField(
        "Additional columns", min=0, default=0, step=1
    )
    structure.general.column.lb1 = LineBreak()
    structure.general.column.profile = OptionField(
        "Column profile", options=profiles_options, default="SHS 50x50 x 4"
    )

    structure.general.purlin = Section("Warehouse roof beams")
    structure.general.purlin.text_ = Text(
        "Purlins are used to distribute the roof load over the trusses."
    )
    structure.general.purlin.profile = OptionField(
        "Purlin profile", options=profiles_options, default="SHS 50x50 x 4"
    )
    structure.general.purlin.lb1 = LineBreak()
    structure.general.purlin.purlin_spacing = OptionField(
        "Purlin per panel", options=PURLIN_SPACING_OPTIONS, default=1.0
    )

    structure.advanced_settings = Tab("Advanced settings")
    structure.advanced_settings.visual_settings = Section("Visibility")
    structure.advanced_settings.visual_settings.office_visible = ToggleButton(
        "Office visible", default=True
    )
    structure.advanced_settings.visual_settings.warehouse_visible = ToggleButton(
        "Warehouse visible", default=True
    )

    final_step = Step("What's next?", views="final_step")
