/*
  Copyright (C) 2026 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#include "catch.hpp"
#include "common.h"
#include <aspect/material_model/utilities.h>
#include <aspect/utilities.h>

namespace
{
  std::string
  lookup_test_file(const std::string &filename)
  {
    return aspect::Utilities::expand_ASPECT_SOURCE_DIR(
             "$ASPECT_SOURCE_DIR/data/material-model/steinberger/test_lookup/"
             + filename);
  }
}

TEST_CASE("Material lookup optional thermal conductivity column")
{
  using aspect::MaterialModel::MaterialUtilities::Lookup::PerplexReader;

  const PerplexReader conductivity_table(lookup_test_file("conductivity.dat"),
                                        true,
                                        MPI_COMM_WORLD);
  REQUIRE(conductivity_table.has_thermal_conductivity());
  CHECK(conductivity_table.thermal_conductivity(350.0, 5e8) == Approx(5.0));
  CHECK(conductivity_table.thermal_conductivity(300.0, 0.0) == Approx(2.0));
  CHECK(conductivity_table.thermal_conductivity(400.0, 1e9) == Approx(8.0));
  CHECK(conductivity_table.thermal_conductivity(350.0, 0.0) == Approx(3.0));
  CHECK(conductivity_table.thermal_conductivity(300.0, 5e8) == Approx(4.0));
  CHECK(conductivity_table.thermal_conductivity(200.0, 0.0) == Approx(2.0));
  CHECK(conductivity_table.thermal_conductivity(500.0, 0.0) == Approx(4.0));
  CHECK(conductivity_table.thermal_conductivity(300.0, -1e8) == Approx(2.0));
  CHECK(conductivity_table.thermal_conductivity(300.0, 2e9) == Approx(6.0));
  CHECK(conductivity_table.seismic_Vp(350.0, 5e8) == Approx(8.0));

  const PerplexReader legacy_table(lookup_test_file("legacy.dat"),
                                   true,
                                   MPI_COMM_WORLD);
  CHECK(legacy_table.has_thermal_conductivity());
  CHECK(legacy_table.thermal_conductivity(350.0, 5e8) == Approx(3.0));

  const PerplexReader legacy_table_without_conductivity(lookup_test_file("legacy-no-conductivity.dat"),
                                                        true,
                                                        MPI_COMM_WORLD);
  CHECK_FALSE(legacy_table_without_conductivity.has_thermal_conductivity());
  CHECK(legacy_table_without_conductivity.seismic_Vp(350.0, 5e8) == Approx(8.0));
  CHECK_THROWS(legacy_table_without_conductivity.thermal_conductivity(350.0, 5e8));
}

TEST_CASE("Material lookup rejects invalid conductivity values")
{
  using aspect::MaterialModel::MaterialUtilities::Lookup::PerplexReader;

  CHECK_THROWS_WITH(PerplexReader(lookup_test_file("conductivity-zero.dat"),
                                  true,
                                  MPI_COMM_WORLD),
                    Catch::Matchers::Contains("finite and strictly positive"));
  CHECK_THROWS_WITH(PerplexReader(lookup_test_file("conductivity-negative.dat"),
                                  true,
                                  MPI_COMM_WORLD),
                    Catch::Matchers::Contains("finite and strictly positive"));
  CHECK_THROWS(PerplexReader(lookup_test_file("conductivity-nan.dat"),
                             true,
                             MPI_COMM_WORLD));
}

TEST_CASE("Material lookup interpolation on larger grid")
{
  using aspect::MaterialModel::MaterialUtilities::Lookup::PerplexReader;

  const PerplexReader table(lookup_test_file("conductivity-3x3.dat"),
                            true,
                            MPI_COMM_WORLD);
  CHECK(table.thermal_conductivity(350.0, 2.5e8) == Approx(4.0));
  CHECK(table.thermal_conductivity(450.0, 7.5e8) == Approx(8.0));
}

TEST_CASE("Material lookup parses conductivity after non-finite table values")
{
  using aspect::MaterialModel::MaterialUtilities::Lookup::PerplexReader;

  const std::string filename = aspect::Utilities::expand_ASPECT_SOURCE_DIR(
                                 "$ASPECT_SOURCE_DIR/data/material-model/steinberger/"
                                 "pyr_MS95_with_volume_fractions_lo_res_k.dat");
  const PerplexReader table(filename, true, MPI_COMM_WORLD);

  CHECK(table.thermal_conductivity(3800.0, 1e9) == Approx(1.394966755));
}

TEST_CASE("Geometric conductivity volume-fraction mixing")
{
  const std::vector<double> fractions = {0.5, 0.5};
  const std::vector<double> conductivities = {4.0, 16.0};
  const auto average = [](const std::vector<double> &volume_fractions,
                          const std::vector<double> &values)
  {
    return aspect::MaterialModel::MaterialUtilities::average_value(
             volume_fractions,
             values,
             aspect::MaterialModel::MaterialUtilities::geometric);
  };

  CHECK(average(fractions, conductivities) == Approx(8.0));
  CHECK_THROWS_WITH(average({1.0}, conductivities),
                    Contains("have to have the same length"));
}
