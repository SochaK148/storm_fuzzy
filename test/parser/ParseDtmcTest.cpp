/*
 * ParseDtmcTest.cpp
 *
 *  Created on: 03.12.2012
 *      Author: Thomas Heinemann
 */


#include "gtest/gtest.h"
#include "storm-config.h"
#include "src/parser/DtmcParser.h"
#include "src/utility/IoUtility.h"

TEST(ParseDtmcTest, parseAndOutput) {
	storm::parser::DtmcParser* dtmcParser;
	ASSERT_NO_THROW(dtmcParser = new storm::parser::DtmcParser(
			STORM_CPP_TESTS_BASE_PATH "/parser/tra_files/pctl_general_input_01.tra",
			STORM_CPP_TESTS_BASE_PATH "/parser/lab_files/pctl_general_input_01.lab"));

	ASSERT_NO_THROW(storm::utility::dtmcToDot(*(dtmcParser->getDtmc()), STORM_CPP_TESTS_BASE_PATH "/parser/output.dot"));

	delete dtmcParser;
}


