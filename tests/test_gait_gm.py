#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `gait_gm` package."""


import unittest
from click.testing import CliRunner

from gait_gm import gait_gm
from gait_gm import cli


class TestGait_gm(unittest.TestCase):
    """Tests for `gait_gm` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_000_something(self):
        """Test something."""

    def test_command_line_interface(self):
        """Test the CLI."""
        runner = CliRunner()
        result = runner.invoke(cli.main)
        assert result.exit_code == 0
        assert 'gait_gm.cli.main' in result.output
        help_result = runner.invoke(cli.main, ['--help'])
        assert help_result.exit_code == 0
        assert '--help  Show this message and exit.' in help_result.output
