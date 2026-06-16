import unittest
import numpy as np
import astropy.units as u
import radiation_unit_conversion.units as units

test_flux = np.linspace(1, 10, 10)
test_wavelength = np.linspace(1, 5, 10)


class TestUnitsConversions(unittest.TestCase):


    def setUp(self):
        """Set up test data"""
        self.test_flux = np.linspace(1, 10, 10)  # Array of flux values from 1 to 10 W/m^2
        self.test_wavelength = np.linspace(1, 5, 10)  # Array of wavelength values from 1 to 5 micron

    def test_watt_metersquared2erg_cmsquared_second(self):
        """Test conversion from W/m^2 to erg/cm^2/s"""
        expected_output = self.test_flux * 1000  # erg/cm^2/s
        result = units.watt_metersquared2erg_cmsquared_second(self.test_flux)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_erg_cmsquared_second2watt_metersquared(self):
        """Test conversion from erg/cm^2/s to W/m^2"""
        test_flux_erg = self.test_flux * 1000  # erg/cm^2/s
        expected_output = self.test_flux  # W/m^2
        result = units.erg_cmsquared_second2watt_metersquared(test_flux_erg)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_round_trip_conversion(self):
        """Test that a round-trip conversion returns the original values."""
        intermediate_flux = units.watt_metersquared2erg_cmsquared_second(self.test_flux)
        final_flux = units.erg_cmsquared_second2watt_metersquared(intermediate_flux)
        np.testing.assert_allclose(self.test_flux, final_flux, rtol=1e-6)

    def test_watt_metersquared_hertz2erg_cmsquared_second_hertz(self):
        """Test conversion from W/m^2/Hz to erg/cm^2/s/Hz"""
        expected_output = self.test_flux * 1000
        result = units.watt_metersquared_hertz2erg_cmsquared_second_hertz(self.test_flux)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_erg_cmsquared_secondhertz2watt_metersquaredhertz(self):
        """Test conversion from erg/cm^2/s/Hz to W/m^2/Hz"""
        test_flux_erg = self.test_flux * 1000
        expected_output = self.test_flux
        result = units.erg_cmsquared_second_hertz2watt_metersquaredhertz(test_flux_erg)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_watt_metersquared_hertz2erg_cmsquared_second_angstrom(self):
        """Test conversion from W/m^2/Hz to erg/cm^2/s/\u00c5ngström"""
        constant = 2.99792458e21
        expected_output = constant * self.test_flux / self.test_wavelength**2
        result = units.watt_metersquared_hertz2erg_cmsquared_second_angstrom(self.test_flux, self.test_wavelength)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_erg_cmsquared_second_angstrom2watt_metersquared_hertz(self):
        """Test conversion from erg/cm^2/s/\u00c5ngström to W/m^2/Hz"""
        test_flux_erg = self.test_flux * 1000
        constant = 2.99792458e21
        expected_output =  test_flux_erg * self.test_wavelength**2 / constant
        result = units.erg_cmsquared_second_angstrom2watt_metersquared_hertz(test_flux_erg, self.test_wavelength)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_watt_metersquared_micron2watt_metersquared_hertz(self):
        """Test conversion from W/m^2/\u03bcm to W/m^2/Hz"""
        constant = 2.99792458e14
        expected_output = self.test_flux * self.test_wavelength**2 / constant
        result = units.watt_metersquared_micron2watt_metersquared_hertz(self.test_flux, self.test_wavelength)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_watt_metersquared_hertz2watt_metersquared_micron(self):
        """Test conversion from W/m^2/Hz to W/m^2/\u03bcm"""
        constant = 2.99792458e14
        test_flux_hz =  self.test_flux * self.test_wavelength**2 / constant
        expected_output = self.test_flux
        result = units.watt_metersquared_hertz2watt_metersquared_micron(test_flux_hz, self.test_wavelength)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_photon_cmsquared_second_micron2watt_metersquared_micron(self):
        """Test conversion from photon/cm^2/s/μm to W/m^2/μm"""
        constant = 5.03411250e14
        expected_output = self.test_flux / self.test_wavelength / constant
        result = units.photon_cmsquared_second_micron2watt_metersquared_micron(self.test_flux, self.test_wavelength)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_erg_cmsquared_second_angstrom2photon_cmsquared_second_angstrom(self):
        """Test conversion from erg/cm^2/s/Ångström to photon/cm^2/s/Ångström"""
        constant = 5.03411250E+07
        expected_output = constant * self.test_flux * self.test_wavelength
        result = units.erg_cmsquared_second_angstrom2photon_cmsquared_second_angstrom(self.test_flux, self.test_wavelength)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_photon_cmsquared_second_angstrom2erg_cmsquared_second_angstrom(self):
        """Test conversion from photon/cm^2/s/Ångström to erg/cm^2/s/Ångström"""
        constant = 5.03411250E+07
        expected_output = self.test_flux / self.test_wavelength / constant
        result = units.photon_cmsquared_second_angstrom2erg_cmsquared_second_angstrom(self.test_flux, self.test_wavelength)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_watt_metersquared_hertz2jansky(self):
        """Test conversion from W/m^2/Hz to Jansky"""
        constant = 1e26
        expected_output = self.test_flux * constant
        result = units.watt_metersquared_hertz2jansky(self.test_flux)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_jansky2watt_metersquared_hertz(self):
        """Test conversion from Jansky to W/m^2/Hz"""
        constant = 1e26
        test_flux_jansky = self.test_flux * constant
        expected_output = self.test_flux
        result = units.jansky2watt_metersquared_hertz(test_flux_jansky)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_erg_cmsquared_second_hertz2jansky(self):
        """Test conversion from erg/cm^2/s/Hz to Jansky"""
        constant = 1e23
        expected_output = self.test_flux * constant
        result = units.erg_cmsquared_second_hertz2jansky(self.test_flux)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_jansky2erg_cmsquared_second_hertz(self):
        """Test conversion from Jansky to erg/cm^2/s/Hz"""
        constant = 1e23
        test_flux_jansky = self.test_flux * constant
        expected_output = self.test_flux
        result = units.jansky2erg_cmsquared_second_hertz(test_flux_jansky)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_erg_cmsquared_second_angstrom2jansky(self):
        """Test conversion from erg/cm^2/s/Ångström to Jansky"""
        constant = 3.33564095e04
        expected_output = constant * self.test_flux * self.test_wavelength**2
        result = units.erg_cmsquared_second_angstrom2jansky(self.test_flux, self.test_wavelength)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_jansky2erg_cmsquared_second_angstrom(self):
        """Test conversion from Jansky to erg/cm^2/s/Ångström"""
        constant = 3.33564095e04
        test_flux_jansky = constant * self.test_flux * self.test_wavelength**2
        expected_output = self.test_flux
        result = units.jansky2erg_cmsquared_second_angstrom(test_flux_jansky, self.test_wavelength)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_watt_metersquared_micron2jansky(self):
        """Test conversion from W/m^2/μm to Jansky"""
        # [Y Jy] = 3.33564095E+11 * [X1 W/m^2/um] * [X2 um]^2
        constant = 3.33564095e11
        expected_output = constant * self.test_flux * self.test_wavelength**2
        result = units.watt_metersquared_micron2jansky(self.test_flux, self.test_wavelength)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_jansky2watt_metersquared_micron(self):
        """Test that Jy -> W/m^2/μm inverts W/m^2/μm -> Jy (round trip)."""
        jansky = units.watt_metersquared_micron2jansky(self.test_flux, self.test_wavelength)
        result = units.jansky2watt_metersquared_micron(jansky, self.test_wavelength)
        np.testing.assert_allclose(result, self.test_flux, rtol=1e-6)

    def test_flambda2fnu(self):
        """Test conversion from W/m^2/μm to Jansky"""
        # f_nu = f_lambda * lambda^2 / c; in Jy this is *1e26/c[um/s] = *3.33564095e11
        constant = 1e26 / 2.99792458e14
        expected_output = self.test_flux * self.test_wavelength**2 * constant
        result = units.flambda2fnu(self.test_flux * u.W / u.m**2 / u.micron,
                                   self.test_wavelength * u.micron,
                                   u.Jy)
        np.testing.assert_allclose(result.value, expected_output, rtol=1e-6)

    def test_fnu2flambda(self):
        """Test conversion from Jansky to W/m^2/μm """
        constant = 1e26 / 2.99792458e14
        fnu = self.test_flux * self.test_wavelength**2 * constant
        result = units.fnu2flambda(fnu * u.Jy,
                                   self.test_wavelength * u.micron,
                                   u.W / u.m**2 / u.micron)
        np.testing.assert_allclose(result.value, self.test_flux, rtol=1e-6)


    def test_rayleigh_arcsec_round_trip(self):
        """Rayleigh -> photon/cm^2/s/A/arcsec^2 -> Rayleigh recovers the input."""
        intermediate = units.rayleigh2photon_cmsquared_second_angstrom_arcsecondsquared(self.test_flux)
        result = units.photon_cmsquared_second_angstrom_arcsecondsquared2rayleigh(intermediate)
        np.testing.assert_allclose(result, self.test_flux, rtol=1e-6)


class TestConversionDispatcher(unittest.TestCase):
    """Exercise the high-level conversion() dispatcher with both astropy and float inputs."""

    def setUp(self):
        self.flux = np.linspace(1, 10, 10)
        self.wavelength = np.linspace(1, 5, 10)

    def test_astropy_no_wavelength(self):
        """Astropy flux with no wavelength should not require an explicit unit string."""
        result = units.conversion(self.flux * u.W / u.m**2, output_units='erg/cm2/s')
        np.testing.assert_allclose(result.value, self.flux * 1000, rtol=1e-6)
        self.assertEqual(result.unit, u.erg / u.cm**2 / u.s)

    def test_float_string_units(self):
        """Float flux with explicit string units."""
        result = units.conversion(self.flux, output_units='erg/cm2/s', in_unit='W/m2')
        np.testing.assert_allclose(result, self.flux * 1000, rtol=1e-6)

    def test_jy_to_flambda_round_trip_astropy(self):
        """W/m^2/um -> Jy -> W/m^2/um using astropy quantities recovers the input."""
        flambda = self.flux * u.W / u.m**2 / u.micron
        wavelength = self.wavelength * u.micron
        jy = units.conversion(flambda, output_units='Jy', in_wavelength=wavelength)
        back = units.conversion(jy, output_units='W/m2/um', in_wavelength=wavelength)
        np.testing.assert_allclose(back.value, self.flux, rtol=1e-6)

    def test_cross_unit_wavelength_float(self):
        """Wavelength supplied in microns for an Angstrom-based conversion (the old NameError path)."""
        result = units.conversion(
            self.flux,
            output_units='photon/cm2/s/A',
            in_unit='erg/cm2/s/A',
            in_wavelength=self.wavelength,
            in_wavelength_units='um',
        )
        wavelength_angstrom = self.wavelength * 1e4
        expected = 5.03411250e07 * self.flux * wavelength_angstrom
        np.testing.assert_allclose(result, expected, rtol=1e-6)

    def test_missing_units_raises(self):
        with self.assertRaises(ValueError):
            units.conversion(self.flux, output_units='erg/cm2/s')

    def test_unsupported_conversion_raises(self):
        with self.assertRaises(ValueError):
            units.conversion(self.flux, output_units='not/a/unit', in_unit='W/m2')


if __name__ == "__main__":
    unittest.main()
