from pathlib import Path

from munch import Munch
from viktor import ViktorController
from viktor.geometry import CartesianAxes, Group, Point
from viktor.views import (
    GeometryResult,
    GeometryView,
    MapResult,
    MapView,
    WebResult,
    WebView,
)

from app.parametrization import WarehouseParametrization

from .model import BuildingExterior, Map, OfficeFrame, WarehouseSteelFrame


class WarehouseController(ViktorController):
    label = "Warehouse configurator"
    parametrization = WarehouseParametrization
    viktor_enforce_field_constraints = True

    @MapView("Map View", duration_guess=1)
    def get_map_view(self, params: Munch, **kwargs):
        """Initiates the process of rendering a Map of the terrain and location of the warehouse model."""
        map_plot = Map.from_params(params=params)

        features = []

        if map_plot.land_polygon:
            features.append(map_plot.get_land_polygon())

        if map_plot.building_corner:
            features.append(map_plot.get_office_polygon())
            features.append(map_plot.get_warehouse_polygon())

        return MapResult(features)

    @GeometryView("Building", duration_guess=1)
    def visualize_exterior(self, params, **kwargs):
        """Initiates the process of rendering a 3D model of the exteriorof the warehouse model."""

        exterior = BuildingExterior.from_params(params).visualize()
        land = Map.from_params(params=params).visualize()

        cartesian = CartesianAxes(Point(0, 0, 0), 10, 0.1)
        exterior = Group([exterior, land, cartesian])

        return GeometryResult(exterior)

    @GeometryView("Structure", duration_guess=10)
    def visualize_structure(self, params, **kwargs):
        """Initiates the process of rendering a the steel structure of the warehouse model."""

        if params.structure.advanced_settings.visual_settings.warehouse_visible:
            warehouse = WarehouseSteelFrame.from_params(params=params).visualise()
        else:
            warehouse = WarehouseSteelFrame.from_params(params).draw_floor()

        if params.structure.advanced_settings.visual_settings.office_visible:
            office = OfficeFrame.from_params(params=params).visualise()
        else:
            office = OfficeFrame.from_params(params=params).draw_floor()

        beams = Group([warehouse, office])

        return GeometryResult(beams)

    @WebView(" ", duration_guess=1)
    def final_step(self, params, **kwargs):
        """Initiates the process of rendering the last step."""
        html_path = Path(__file__).parent / "final_step.html"
        with html_path.open() as f:
            html_string = f.read()
        return WebResult(html=html_string)
