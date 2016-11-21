/*
 * SUMMARY:      settings.h - Definition of string, array sizes, etc.
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * LAST-MOD:     02/11/2013 by Ning Sun
 * DESCRIPTION:  Definition of string, array sizes, etc.
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:     
 * $Id: settings.h,v 3.1.1 2013/02/11 ning Exp $	 */

#ifndef SETTINGS_H
#define SETTINGS_H

#ifndef TRUE
#define TRUE           1
#endif
#ifndef FALSE
#define FALSE          0
#endif

#define DHSVM_HUGE     (1e20)

/* When compiling DHSVM using Microsoft C++ define MSC++ at compile time.
   Microsoft C++ treats all numeric constants as doubles, and issues a
   warning if this constant is assigned to a float without an explicit type
   cast.   This warning can becme somewhat annoying, and the following pragma
   directive disables the warning */
#ifdef MS_C_PLUSPLUS
#pragma warning(disable : 4244)
#endif

/* Default value for not applicable */
#define NOT_APPLICABLE -9999

/* Options for precipitation and wind source */
#define RADAR          1
#define STATION        2
#define MODEL          3

/* Options for flow gradient calculation */
#define TOPOGRAPHY     1
#define WATERTABLE     2

/* Options for model extent */
#define POINT 1
#define BASIN 2

/* Options for temperature and precipitation lapse rates */
#define CONSTANT 1
#define VARIABLE 2

/* Options for radiation */
#define IPW     1               /* use IPW to calculate solar radiation */
#define INLINE  2               /* use inline calculations to calculate 
                                      solar radiation for each pixel */ 
 
#define TINY       1e-20
#define DEBUG      FALSE

#define INFOFILE "input.filenames.txt"

#define HEADERLINES    5
#define BUFSIZE      255
#define MAXUCHAR     255	/* Maximum value of a 1-byte uchar */
#define MAXSTRING    255
#define NAMESIZE     127

#define NDIRS          4	/* Number of directions in which water can 
				   flow */
#define NA          -999	/* Not applicable */

#define N_MM5_MAPS 7

#define MAP_OUTPUT 1
#define IMAGE_OUTPUT 2

enum KEYS {
  /* Options */
  format = 0, extent, gradient, flow_routing, sensible_heat_flux, 
  mm5, radiation, precipitation_source, wind_source, temp_lapse,
  precip_lapse,
  /* Area */
  coordinate_system, extreme_north, extreme_west, center_latitude,
  center_longitude, time_zone_meridian, number_of_rows, 
  number_of_columns, grid_spacing, point_north, point_east, 
  /* Time */
  time_step, model_start, model_end, 
  /* Constants */
  ground_roughness, snow_roughness, rain_threshold, snow_threshold, 
  snow_water_capacity, reference_height, rain_lai_multiplier, 
  snow_lai_multiplier, min_intercepted_snow, outside_basin, 
  temp_lapse_rate, precip_lapse_rate,
  /* Station information */
  station_name = 0, station_north, station_east, station_elev, station_file,
  /* IPW information */
  IPW_table_file = 0, IPW_map_file,
  /* RADAR information */
  radar_start = 0, radar_file, radar_north, radar_west, radar_rows, radar_cols,
  radar_grid,
  /* Wind model information */
  number_of_maps = 0, wind_map_path, wind_station,
  /* MM5 information */
  MM5_start = 0, MM5_temperature, MM5_humidity, MM5_wind, MM5_shortwave,
  MM5_longwave, MM5_pressure, MM5_precip,
  /* Soil information */
  soil_description = 0, lateral_ks, exponent, infiltration, soil_albedo, 
  number_of_layers, porosity, pore_size, bubbling_pressure, field_capacity, 
  wilting_point, bulk_density, vertical_ks, solids_thermal, thermal_capacity,
  /* Vegetation information */
  veg_description = 0, overstory, understory, fraction, trunk_space, 
  aerodynamic_att, radiation_att, height, summer_lai, winter_lai,
  max_resistance, min_resistance, moisture_threshold, vpd, rpc,
  veg_albedo, number_of_root_zones, root_zone_depth, overstory_fraction, 
  understory_fraction,
  /* terrain information */
  demfile = 0, maskfile, 
  soiltype_file = 0, soildepth_file,
  /* DHSVM channel keys */
  stream_network = 0, stream_map, stream_class,  
  road_network, road_map, road_class, 
  /* number of each type of output */
  output_path = 0, initial_state_path, npixels, nstates, nmapvars, nimagevars,
  /* pixel information */
  north = 0, east,
  /* state information */
  state_date = 0,
  /* map information */
  map_variable = 0, map_layer, nmaps, map_date,
  /* image information */
  image_variable = 0, image_layer, image_start, image_end, image_interval,
  image_upper, image_lower
};

#endif

