from dataclasses import dataclass
from math import ceil, cos, pi, sin
from typing import Tuple

import geopy.distance
import numpy as np
import viktor as vkt

from .database import profile_properties

FLOOR_HEIGHT = 4  # [m]
WINDOWS_HEIGHT = 2  # [m]
FACADE_HEIGHT = FLOOR_HEIGHT - WINDOWS_HEIGHT  # [m]


class BuildingExterior:
    def __init__(
        self,
        width_building,
        width_office,
        num_office_floors,
        width_warehouse,
        num_warehouse_floors,
        **kwargs
    ):
        # Office properties
        self.width_building = width_building
        self.width_office = width_office
        self.num_office_floors = num_office_floors

        # Warehouse properties
        self.width_warehouse = width_warehouse
        self.num_warehouse_floors = num_warehouse_floors

    @property
    def facade_height(self):
        return FLOOR_HEIGHT - WINDOWS_HEIGHT

    @property
    def height_office(self):
        return self.num_office_floors * FLOOR_HEIGHT + self.facade_height

    @property
    def height_warehouse(self):
        return self.num_warehouse_floors * FLOOR_HEIGHT + self.facade_height / 2

    @property
    def facade_material(self):
        return vkt.Material("Facade", color=vkt.Color(200, 200, 255), opacity=1, metalness=0.5)

    @property
    def windows_material(self):
        return vkt.Material("Windows", color=vkt.Color(180, 180, 210), opacity=0.90, metalness=1)

    @property
    def roof_material(self):
        return vkt.Material("Facade", color=vkt.Color(74, 80, 84), opacity=1, metalness=0.5)

    @classmethod
    def from_params(cls, params):
        return cls(
            width_building=params.location_and_building.building_y_dimension,
            width_office=params.location_and_building.office_x_dimension,
            num_office_floors=params.location_and_building.num_office_floors,
            width_warehouse=params.location_and_building.warehouse_x_dimension,
            num_warehouse_floors=params.location_and_building.num_warehouse_floors,
        )

    def visualize(self):
        """Returns an object of the building that can be rendered in the `GeometryView`."""
        office = vkt.Group(
            [
                self.draw_office_first_floor(),
                self.draw_office_facades(),
                self.draw_office_windows(),
                self.draw_office_roof(),
            ]
        )

        office.translate([self.width_office / 2, self.width_building / 2, 0])

        warehouse = vkt.Group(
            [
                self.draw_warehouse_facades(),
                self.draw_warehouse_windows(),
                self.draw_warehouse_roof(),
            ]
        )

        # translation done to put building in centre of view.
        warehouse.translate(
            [self.width_office + self.width_warehouse / 2, self.width_building / 2, 0]
        )

        return vkt.Group([office, warehouse])

    def draw_office_facades(self) -> vkt.Group:
        loc_z = FLOOR_HEIGHT
        line = vkt.Line(vkt.Point(0, 0, loc_z), vkt.Point(0, 0, loc_z + self.facade_height))
        facade = vkt.RectangularExtrusion(self.width_office, self.width_building, line)
        facade.material = self.facade_material
        pattern = vkt.LinearPattern(facade, [0, 0, 1], self.num_office_floors, FLOOR_HEIGHT)
        return pattern

    def draw_office_windows(self) -> vkt.Group:
        loc_z = FLOOR_HEIGHT + self.facade_height
        line = vkt.Line(vkt.Point(0, 0, loc_z), vkt.Point(0, 0, loc_z + WINDOWS_HEIGHT))
        window = vkt.RectangularExtrusion(self.width_office, self.width_building, line)
        window.material = self.windows_material
        pattern = vkt.LinearPattern(window, [0, 0, 1], self.num_office_floors - 1, FLOOR_HEIGHT)
        return pattern

    def draw_office_first_floor(self):
        line = vkt.Line(vkt.Point(0, 0, 0), vkt.Point(0, 0, FLOOR_HEIGHT))
        window = vkt.RectangularExtrusion(self.width_office, self.width_building, line)
        window.material = self.windows_material
        return window

    def draw_office_roof(self):
        line = vkt.Line(vkt.Point(0, 0, self.height_office), vkt.Point(0, 0, self.height_office + 0.5))
        roof = vkt.RectangularExtrusion(self.width_office, self.width_building, line)
        roof.material = self.roof_material
        return roof

    def draw_warehouse_facades(self) -> vkt.Group:
        loc_z = self.height_warehouse - WINDOWS_HEIGHT - self.facade_height / 2

        line = vkt.Line(vkt.Point(0, 0, 0), vkt.Point(0, 0, loc_z))
        facade1 = vkt.RectangularExtrusion(self.width_warehouse, self.width_building, line)
        facade1.material = self.facade_material

        line = vkt.Line(
            vkt.Point(0, 0, loc_z + WINDOWS_HEIGHT), vkt.Point(0, 0, self.height_warehouse)
        )
        facade2 = vkt.RectangularExtrusion(self.width_warehouse, self.width_building, line)
        facade2.material = self.facade_material

        return vkt.Group([facade1, facade2])

    def draw_warehouse_windows(self) -> vkt.Group:
        loc_z = self.height_warehouse - WINDOWS_HEIGHT - self.facade_height / 2
        line = vkt.Line(vkt.Point(0, 0, loc_z), vkt.Point(0, 0, loc_z + WINDOWS_HEIGHT))
        window = vkt.RectangularExtrusion(self.width_warehouse, self.width_building, line)
        window.material = self.windows_material
        return window

    def draw_warehouse_roof(self):
        line = vkt.Line(vkt.Point(0, 0, self.height_warehouse), vkt.Point(0, 0, self.height_warehouse + 0.5))
        roof = vkt.RectangularExtrusion(self.width_warehouse, self.width_building, line)
        roof.material = self.roof_material
        return roof


@dataclass
class BeamProfile:
    width: float
    mass: float
    area: float
    inertia: float
    section_modulus: float


class WarehouseSteelFrame:
    def __init__(
        self,
        y_dimension: float,
        x_dimension: float,
        num_warehouse_floors: int,
        num_columns: int,
        column_profile: BeamProfile,
        truss_max_spacing: float,
        truss_height: float,
        truss_chord_profile: BeamProfile,
        truss_web_profile: BeamProfile,
        truss_vertical_profile: BeamProfile,
        purlin_spacing: float,
        purlin_profile: BeamProfile,
        num_truss_panels: int = 0,
        **kwargs
    ):
        self.material = vkt.Material("blue", color=vkt.Color(100, 100, 255))

        # Warehouse properties
        self.width_warehouse = x_dimension  # Frame spacing direction
        self.width_building = y_dimension  # Truss direction
        self.warehouse_height = num_warehouse_floors * FLOOR_HEIGHT + FACADE_HEIGHT / 2

        # Columns properties
        self.num_columns = num_columns + 2
        self.column_profile = column_profile

        # Frames properties
        self.num_frames = ceil(self.width_warehouse / truss_max_spacing) + 1
        self.frame_spacing = self.width_warehouse / (self.num_frames - 1)
        self.num_truss_per_frame = self.num_columns - 1

        # Truss properties
        self.truss_length = self.width_building / self.num_truss_per_frame
        self.truss_height = truss_height

        if num_truss_panels:
            self.truss_panels = 2 * int(num_truss_panels / 2)
        else:
            self.truss_panels = 2 * int(
                self.width_building / (2 * self.truss_height * self.num_truss_per_frame)
            )

        self.truss_chord = truss_chord_profile
        self.truss_web = truss_web_profile
        self.truss_vertical = truss_vertical_profile

        # Purlin
        self.purlin_spacing = purlin_spacing
        self.purlin_profile = purlin_profile

        # Floor height:
        self.size_z = self.warehouse_height - (
            self.truss_height + 0.5 * self.truss_chord.width + self.purlin_profile.width
        )

        self.direction_x = [1, 0, 0]
        self.direction_y = [0, 1, 0]

    @classmethod
    def from_params(cls, params):
        column_profile = BeamProfile(**profile_properties[params.structure.general.column.profile])
        purlin_profile = BeamProfile(**profile_properties[params.structure.general.purlin.profile])
        truss_chord_profile = BeamProfile(**profile_properties[params.structure.general.truss.profile_chord])
        truss_web_profile = BeamProfile(**profile_properties[params.structure.general.truss.profile_web])
        truss_vertical_profile = BeamProfile(**profile_properties[params.structure.general.truss.profile_vertical])
        custom_panels = params.structure.general.truss.custom_panels
        return cls(
            y_dimension=params.location_and_building.building_y_dimension,
            x_dimension=params.location_and_building.warehouse_x_dimension,
            num_warehouse_floors=params.location_and_building.num_warehouse_floors,
            num_columns=params.structure.general.column.num_columns,
            column_profile=column_profile,
            truss_max_spacing=params.structure.general.truss.max_truss_spacing,
            truss_height=params.structure.general.truss.truss_height,
            truss_chord_profile=truss_chord_profile,
            truss_web_profile=truss_web_profile,
            truss_vertical_profile=truss_vertical_profile,
            num_truss_panels=params.structure.general.truss.truss_panels if custom_panels else None,
            purlin_spacing=params.structure.general.purlin.purlin_spacing,
            purlin_profile=purlin_profile,
        )

    def visualise(self):
        return vkt.Group([
            self.draw_trusses(),
            self.draw_columns(),
            self.draw_purlin(),
            self.draw_floor(),
            self.draw_bracing(),
        ])

    def get_truss(self):
        return Truss(
            length=self.truss_length,
            height=self.truss_height,
            sections=self.truss_panels,
            truss_chord=self.truss_chord,
            truss_web=self.truss_web,
            truss_vertical=self.truss_vertical,
            material=self.material,
        )

    def draw_trusses(self):
        truss = self.get_truss().visualise()

        loc_z_truss_bottom = (
            self.warehouse_height - self.truss_height - self.purlin_profile.width / 2
        )

        truss.rotate(angle=90 * np.pi / 180, direction=[0, 0, 1]).translate(
            (0, 0, loc_z_truss_bottom)
        )

        pattern = vkt.BidirectionalPattern(
            base_object=truss,
            direction_1=self.direction_x,
            direction_2=self.direction_y,
            number_of_elements_1=self.num_frames,
            number_of_elements_2=self.num_truss_per_frame,
            spacing_1=self.frame_spacing,
            spacing_2=self.truss_length,
        )

        return pattern

    def draw_columns(self) -> vkt.Group:
        column = self.get_column()

        line = vkt.Line(column.start, column.end)

        column = vkt.RectangularExtrusion(column.width, column.width, line)
        column.material = self.material

        pattern = vkt.BidirectionalPattern(
            column,
            self.direction_x,
            self.direction_y,
            self.num_frames,
            self.num_columns,
            self.frame_spacing,
            self.truss_length,
        )

        return pattern

    def get_column(self):
        start = vkt.Point(0, 0, 0)
        end = vkt.Point(0, 0, self.warehouse_height - self.purlin_profile.width)
        return Beam(start, end, self.column_profile)

    def get_purlin(self):
        size = self.purlin_profile.width
        start = vkt.Point(0, 0, self.warehouse_height - size / 2)
        end = vkt.Point(self.width_warehouse, 0, self.warehouse_height - size / 2)
        return Beam(start, end, self.purlin_profile)

    def draw_purlin(self) -> vkt.Group:
        purlin = self.get_purlin()
        purlin_num = int(
            self.num_truss_per_frame * self.truss_panels / self.purlin_spacing
        )
        purlin_spacing = self.width_building / purlin_num

        line = vkt.Line(purlin.start, purlin.end)
        purlin = vkt.RectangularExtrusion(purlin.width, purlin.width, line)
        purlin.material = self.material
        pattern = vkt.LinearPattern(purlin, self.direction_y, purlin_num + 1, purlin_spacing)

        return pattern

    def draw_bracing(self):
        width = 0.150
        thickness = 0.010

        num_height = round(self.warehouse_height / self.frame_spacing)
        size_height = self.warehouse_height / num_height

        num_roof = round(self.width_building / self.frame_spacing)
        size_roof = self.width_building / num_roof

        # wall bracing

        line = vkt.Line(vkt.Point(0, 0, 0), vkt.Point(self.frame_spacing, 0, size_height))
        diag1 = vkt.RectangularExtrusion(width, thickness, line)
        diag1.material = self.material

        line = vkt.Line(vkt.Point(self.frame_spacing, 0, 0), vkt.Point(0, 0, size_height))
        diag2 = vkt.RectangularExtrusion(width, thickness, line)
        diag2.material = self.material

        bracing = vkt.Group([diag1, diag2])
        bracing.add(
            vkt.geometry.mirror_object(bracing, vkt.Point(self.width_warehouse / 2, 0, 0), (1, 0, 0))
        )
        bracing.add(
            vkt.geometry.mirror_object(bracing,vkt.Point(0, self.width_building / 2, 0), (0, 1, 0))
        )

        if num_height > 1:
            pattern1 = vkt.LinearPattern(bracing, [0, 0, 1], num_height, size_height)
        else:
            pattern1 = bracing

        # roof bracing diagonals
        line = vkt.Line(
            vkt.Point(0, 0, self.warehouse_height),
            vkt.Point(0 + self.frame_spacing, size_roof, self.warehouse_height),
        )
        diag1 = vkt.RectangularExtrusion(thickness, width, line)
        diag1.material = self.material

        line = vkt.Line(
            vkt.Point(0 + self.frame_spacing, 0, self.warehouse_height),
            vkt.Point(0, size_roof, self.warehouse_height),
        )
        diag2 = vkt.RectangularExtrusion(thickness, width, line)
        diag2.material = self.material

        bracing = vkt.Group([diag1, diag2])
        bracing.add(
            vkt.geometry.mirror_object(bracing, vkt.Point(self.width_warehouse / 2, 0, 0), [1, 0, 0])
        )
        pattern2 = vkt.LinearPattern(bracing, self.direction_y, num_roof, size_roof)

        return vkt.Group([pattern1, pattern2])

    def draw_floor(self):
        floor_mat = vkt.Material("Floor", color=vkt.Color(180, 180, 180))

        line = vkt.Line(
            vkt.Point(self.width_warehouse / 2, self.width_building / 2, -0.1),
            vkt.Point(self.width_warehouse / 2, self.width_building / 2, 0.1),
        )
        floor = vkt.RectangularExtrusion(self.width_warehouse, self.width_building, line)
        floor.material = floor_mat
        return floor


class OfficeFrame:
    def __init__(
        self,
        y_dimension,
        x_dimension,
        num_office_floors,
        col_spacing_x,
        col_spacing_y,
        col_profile,
        **kwargs
    ):

        self.material = vkt.Material("blue", color=vkt.Color(100, 100, 255))
        self.col_profile = col_profile

        # Office properties
        self.y_dimension = y_dimension
        self.x_dimension = x_dimension
        self.num_office_floors = num_office_floors

        self.facade_height = FLOOR_HEIGHT - WINDOWS_HEIGHT

        self.height = self.num_office_floors * FLOOR_HEIGHT

        self.size = self.col_profile.width

        # pattern grid
        self.col_spacing_y = col_spacing_y
        self.col_spacing_x = col_spacing_x

        self.direction_L = [0, 1, 0]
        self.direction_W = [1, 0, 0]
        self.direction_H = [0, 0, 1]

    @property
    def column_num_x(self):
        return int(self.x_dimension / self.col_spacing_x) + 1

    @property
    def column_num_y(self):
        return int(self.y_dimension / self.col_spacing_y) + 1

    @property
    def grid_size_x_direction(self):
        return self.x_dimension / int(self.x_dimension / self.col_spacing_x)

    @property
    def grid_size_y_direction(self):
        return self.y_dimension / int(self.y_dimension / self.col_spacing_y)

    @classmethod
    def from_params(cls, params):
        col_profile = BeamProfile(**profile_properties[params.structure.general.office_col.col_profile])
        return cls(
            y_dimension=params.location_and_building.building_y_dimension,
            x_dimension=params.location_and_building.office_x_dimension,
            num_office_floors=params.location_and_building.num_office_floors,
            col_spacing_x=params.structure.general.office_col.dist_width,
            col_spacing_y=params.structure.general.office_col.dist_length,
            col_profile=col_profile,
        )

    def visualise(self):
        """Renders the office structure for the `GeometryView`."""
        beams_x = self.draw_beams_x_orientation()
        beams_y = self.draw_beams_y_orientation()
        columns = self.draw_office_columns()
        bracing = self.draw_bracing()
        floor = self.draw_floor()

        return vkt.Group([beams_x, beams_y, columns, bracing, floor])

    def draw_office_columns(self):
        """Generates the vertical column pattern for the office."""
        line = vkt.Line(vkt.Point(-self.x_dimension, 0, 0), vkt.Point(-self.x_dimension, 0, self.height))
        column = vkt.RectangularExtrusion(self.size, self.size, line)
        column.material = self.material
        pattern = vkt.BidirectionalPattern(
            column,
            self.direction_W,
            self.direction_L,
            self.column_num_x,
            self.column_num_y,
            self.grid_size_x_direction,
            self.grid_size_y_direction,
        )
        return pattern

    def draw_beams_y_orientation(self):
        """Generates the horizontal beam pattern with y-direction orientation for the office."""
        line = vkt.Line(
            vkt.Point(-self.x_dimension, 0, FLOOR_HEIGHT),
            vkt.Point(-self.x_dimension, self.y_dimension, FLOOR_HEIGHT),
        )
        column = vkt.RectangularExtrusion(self.size, self.size, line)
        column.material = self.material

        pattern = vkt.BidirectionalPattern(
            column,
            self.direction_W,
            self.direction_H,
            self.column_num_x,
            self.num_office_floors,
            self.grid_size_x_direction,
            FLOOR_HEIGHT,
        )

        return pattern

    def draw_beams_x_orientation(self):
        """Generates the horizontal beam pattern with x-direction orientation for the office."""
        line = vkt.Line(vkt.Point(0, 0, FLOOR_HEIGHT), vkt.Point(-self.x_dimension, 0, FLOOR_HEIGHT))
        column = vkt.RectangularExtrusion(self.size, self.size, line)
        column.material = self.material

        pattern = vkt.BidirectionalPattern(
            column,
            self.direction_L,
            self.direction_H,
            self.column_num_y,
            self.num_office_floors,
            self.grid_size_y_direction,
            FLOOR_HEIGHT,
        )

        return pattern

    def draw_bracing(self):
        """Generates the cross braces pattern for the office."""
        width = 0.150
        thickness = 0.010

        # Front diagonals
        line = vkt.Line(
            vkt.Point(-self.x_dimension, 0, 0),
            vkt.Point(-self.x_dimension, self.grid_size_y_direction, FLOOR_HEIGHT),
        )
        diag1 = vkt.RectangularExtrusion(thickness, width, line)
        diag1.material = self.material

        line = vkt.Line(
            vkt.Point(-self.x_dimension, self.grid_size_y_direction, 0),
            vkt.Point(-self.x_dimension, 0, FLOOR_HEIGHT),
        )
        diag2 = vkt.RectangularExtrusion(thickness, width, line)
        diag2.material = self.material

        # Side diagonals
        line = vkt.Line(vkt.Point(0, 0, 0), vkt.Point(-self.grid_size_x_direction, 0, FLOOR_HEIGHT))
        diag3 = vkt.RectangularExtrusion(width, thickness, line)
        diag3.material = self.material

        line = vkt.Line(vkt.Point(-self.grid_size_x_direction, 0, 0), vkt.Point(0, 0, FLOOR_HEIGHT))
        diag4 = vkt.RectangularExtrusion(width, thickness, line)
        diag4.material = self.material

        diagonals = vkt.Group([diag1, diag2, diag3, diag4])
        diagonals.add(
            vkt.geometry.mirror_object(diagonals, vkt.Point(0, self.y_dimension / 2, 0), (0, 1, 0))
        )

        pattern = vkt.LinearPattern(
            diagonals, self.direction_H, self.num_office_floors, FLOOR_HEIGHT
        )

        return pattern

    def draw_floor(self):
        """Generates the floor of the office."""
        floor_mat = vkt.Material("Floor", color=vkt.Color(180, 180, 180))

        line = vkt.Line(
            vkt.Point(-self.x_dimension / 2, self.y_dimension / 2, -0.1),
            vkt.Point(-self.x_dimension / 2, self.y_dimension / 2, 0.1),
        )
        floor = vkt.RectangularExtrusion(self.x_dimension, self.y_dimension, line)
        floor.material = floor_mat

        return floor


class Map:
    def __init__(
        self,
        land_polygon: vkt.GeoPolygon,
        building_corner_coordinate: vkt.GeoPoint,
        building_rotation: float,
        building_y_dimension: float,
        office_x_dimension: float,
        warehouse_x_dimension: float,
        **kwargs
    ):
        # Land
        self.land_polygon = land_polygon

        # Building location
        self.building_corner = building_corner_coordinate
        self.building_rotation = building_rotation * pi / 180

        # Building dimensions
        self.building_y_dimension = (
            self.office_y_dimension
        ) = self.warehouse_y_dimension = building_y_dimension
        self.office_x_dimension = office_x_dimension
        self.warehouse_x_dimension = warehouse_x_dimension

    @classmethod
    def from_params(cls, params):
        return cls(
            land_polygon=params.location_and_building.poly,
            building_corner_coordinate=params.location_and_building.start,
            building_rotation=params.location_and_building.rotate,
            building_y_dimension=params.location_and_building.building_y_dimension,
            office_x_dimension=params.location_and_building.office_x_dimension,
            warehouse_x_dimension=params.location_and_building.warehouse_x_dimension,
        )

    def get_land_polygon(self):
        """Returns a `MapPolygon` that can display the region footprint in a `MapView`."""
        return vkt.MapPolygon.from_geo_polygon(self.land_polygon)

    def visualize(self):
        """Returns an object that displays the selected region in the `GeometryView`."""
        land_points = self._convert_map_coordinates_to_cartesian()
        profile = []

        for point in land_points:
            profile.append(vkt.Point(point[0], point[1], 0))

        profile.append(profile[0])

        if vkt.geometry.circumference_is_clockwise(profile) is not True:
            profile.reverse()

        line = vkt.Line(vkt.Point(0, 0, -0.30), vkt.Point(0, 0, -0.1))
        land = vkt.Extrusion(profile, line).rotate(-self.building_rotation, [0, 0, 1])
        land.material = vkt.Material(
            "grass", color=vkt.Color(120, 255, 130), opacity=0.9
        )

        return land

    def _convert_map_coordinates_to_cartesian(self):
        origin = self.building_corner

        geo_points_list = self.land_polygon.points

        points = []
        for coordinate in geo_points_list:
            x = geopy.distance.distance(
                (origin.lat, origin.lon), (origin.lat, coordinate.lon)
            ).m
            y = geopy.distance.distance(
                (origin.lat, origin.lon), (coordinate.lat, origin.lon)
            ).m

            if coordinate.lon < origin.lon:
                x = -x
            if coordinate.lat < origin.lat:
                y = -y

            points.append((x, y))

        return points

    def get_office_polygon(self):
        """Returns the footprint map polygon of the office part of the building."""
        y_width = self.office_y_dimension
        x_width = self.office_x_dimension
        return self._get_building_polygon(y_width, x_width)

    def get_warehouse_polygon(self):
        """Returns the footprint map polygon of the warehouse part of the building."""
        y_width = self.warehouse_y_dimension
        x_width = self.warehouse_x_dimension
        return self._get_building_polygon(y_width, x_width, translate=(self.office_x_dimension, 0))

    def _get_building_polygon(self, y_dimension: float, x_dimension: float, translate: Tuple = (0, 0)):
        """
        Returns the building footprint as a polygon for a Map visualization.

        :param y_dimension: The longitudinal dimension of the building, assuming zero rotation
        :param x_dimension: The latitudinal dimension of the building, assuming zero rotation
        :param translate: Translational dimensions, before rotation.
        :return: `MapPolygon` footprint of building
        """
        rotation_angle = self.building_rotation
        start = geopy.Point(self.building_corner.lat, self.building_corner.lon)
        shape_points = self.get_shape_coordinates(y_dimension, x_dimension, translate)
        rotated_shape_points = [
            self.rotate((0, 0), pnt, rotation_angle) for pnt in shape_points
        ]
        coordinates = self._convert_points_to_map_coordinates(
            start, rotated_shape_points
        )
        return vkt.MapPolygon(
            points=[vkt.MapPoint(coord[0], coord[1]) for coord in coordinates],
            color=vkt.Color(0, 0, 255),
        )

    @staticmethod
    def rotate(origin, point, angle):
        """
        Rotate a point counterclockwise by a given angle around a given origin.

        The angle should be given in radians.
        """
        ox, oy = origin
        px, py = point

        qx = ox + cos(angle) * (px - ox) - sin(angle) * (py - oy)
        qy = oy + sin(angle) * (px - ox) + cos(angle) * (py - oy)
        return qx, qy

    def _convert_points_to_map_coordinates(self, start: tuple, points: list):
        start_point = geopy.Point(*start)
        map_coordinates = [self._get_point(start_point, points[0])]
        for i, point in enumerate(points):
            if i == 0:
                continue
            _start_point = map_coordinates[i - 1]
            _x_diff = points[i][0] - points[i - 1][0]
            _y_diff = points[i][1] - points[i - 1][1]
            next_point = self._get_point(_start_point, (_x_diff, _y_diff))
            map_coordinates.append(next_point)
        return map_coordinates

    @staticmethod
    def _get_point(point: geopy.Point, translate: Tuple):
        _x_diff, _y_diff = translate
        _distance = np.sqrt(_x_diff**2 + _y_diff**2)
        d = geopy.distance.distance(kilometers=_distance / 1000)
        if _y_diff == 0:
            if _x_diff >= 0:
                bearing = 90
            else:
                bearing = -90
        else:
            bearing = np.rad2deg(np.arctan(_x_diff / _y_diff))
            if _y_diff < 0:
                bearing = 180 + bearing
        return d.destination(point=point, bearing=bearing)

    @staticmethod
    def get_shape_coordinates(y_dimension, x_dimension, translate: Tuple = (0, 0)):
        start_x, start_y = translate

        point1 = (start_x, start_y)
        point2 = (start_x, start_y + y_dimension)
        point3 = (start_x + x_dimension, start_y + y_dimension)
        point4 = (start_x + x_dimension, start_y)
        return [point1, point2, point3, point4]


class Truss:
    def __init__(
        self, length, height, sections, truss_chord, truss_web, truss_vertical, material
    ):

        self.length = length
        self.sections = sections

        self.material = material

        self.truss_chord = truss_chord
        self.truss_web = truss_web
        self.truss_vertical = truss_vertical

        self.chord_width = truss_chord.width
        self.web_width = truss_web.width
        self.vertical_width = truss_vertical.width

        self.height = height - self.chord_width

        # Components
        self.nodes_bottom = self._create_bottom_nodes()
        self.nodes_top = self._create_top_nodes()
        self.top_beams = self._create_top_beams()
        self.bottom_beams = self._create_bottom_beams()
        self.vertical_beams = self._create_vertical_beams()
        self.diagonal_beams = self._create_diagonal_beams()

        self.beams = [
            *self.top_beams,
            *self.bottom_beams,
            *self.vertical_beams,
            *self.diagonal_beams,
        ]

        # Calculated properties:
        self.weight = sum(beam.weight for beam in self.beams)

    def _create_top_nodes(self) -> list:
        nodes_top = []

        for loc in np.linspace(0, self.length, self.sections + 1):
            nodes_top.append(vkt.Point(loc, 0, self.height))

        return nodes_top

    def _create_bottom_nodes(self) -> list:
        nodes_bottom = []

        for loc in np.linspace(0, self.length, self.sections + 1):
            nodes_bottom.append(vkt.Point(loc, 0, 0))

        return nodes_bottom

    def _create_top_beams(self):
        beams_top = []

        for i in range(len(self.nodes_bottom)):
            if i > 0:
                beams_top.append(
                    Beam(self.nodes_top[i - 1], self.nodes_top[i], self.truss_chord)
                )

        return beams_top

    def _create_bottom_beams(self):
        beams_bottom = []

        for i in range(len(self.nodes_bottom)):
            if i > 0:
                beams_bottom.append(
                    Beam(
                        self.nodes_bottom[i - 1], self.nodes_bottom[i], self.truss_chord
                    )
                )

        return beams_bottom

    def _create_vertical_beams(self):
        beams_vertical = []

        for i in range(len(self.nodes_bottom)):

            if i > 0:

                if i < len(self.nodes_bottom) - 1:
                    beams_vertical.append(
                        Beam(
                            self.nodes_bottom[i], self.nodes_top[i], self.truss_vertical
                        )
                    )

        return beams_vertical

    def _create_diagonal_beams(self):
        diagonal_beams = []
        mid_node_index = int(len(self.nodes_bottom) / 2) + 1

        for i in range(len(self.nodes_bottom)):
            if i > 0:
                if i < mid_node_index:
                    diagonal_beams.append(
                        Beam(
                            self.nodes_bottom[i - 1], self.nodes_top[i], self.truss_web
                        )
                    )
                else:
                    diagonal_beams.append(
                        Beam(
                            self.nodes_bottom[i], self.nodes_top[i - 1], self.truss_web
                        )
                    )

        return diagonal_beams

    def visualise(self) -> vkt.Group:
        beam_objects = []

        for beam_top, beam_bottom in zip(self.top_beams, self.bottom_beams):

            beam_objects.append(
                vkt.RectangularExtrusion(
                    self.chord_width,
                    self.chord_width,
                    vkt.Line(beam_top.start, beam_top.end),
                )
            )
            beam_objects.append(
                vkt.RectangularExtrusion(
                    self.chord_width,
                    self.chord_width,
                    vkt.Line(beam_bottom.start, beam_bottom.end),
                )
            )

        for beam in self.diagonal_beams:
            beam_objects.append(
                vkt.RectangularExtrusion(
                    self.web_width, self.web_width, vkt.Line(beam.start, beam.end)
                )
            )

        for beam in self.vertical_beams:
            beam_objects.append(
                vkt.RectangularExtrusion(
                    self.vertical_width, self.vertical_width, vkt.Line(beam.start, beam.end)
                )
            )

        for obj in beam_objects:
            obj.material = self.material

        return vkt.Group(beam_objects)


class Beam:
    def __init__(self, start: vkt.Point, end: vkt.Point, profile: BeamProfile):
        self.start = start
        self.end = end
        self.profile = profile

        # Strength check input
        self.safety_factor = 1.5
        self.yield_stress = 235e6  # Pa
        self.young_modulus = 210e9  # Pa

    @property
    def area(self):
        return self.profile.area

    @property
    def EA(self):
        return self.young_modulus * self.area

    @property
    def EI(self):
        return self.young_modulus * self.inertia

    @property
    def length(self):
        return self.start.vector_to(self.end).magnitude

    @property
    def width(self):
        return self.profile.width

    @property
    def weight(self):
        return self.length * self.profile.mass

    @property
    def inertia(self):
        return self.profile.inertia

    @property
    def section_modulus(self):
        return self.profile.section_modulus
