from munch import Munch

from viktor import ViktorController
from viktor.views import GeometryView, GeometryResult, MapView, MapResult

from app.parametrization import WarehouseParametrization
from viktor.geometry import Group, Point, CartesianAxes
from .model import WarehouseSteelFrame, BuildingExterior, Map, OfficeFrame


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

    @GeometryView("Structure", duration_guess=1)
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
