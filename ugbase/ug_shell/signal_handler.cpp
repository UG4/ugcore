/*
 * signal_handler.cpp
 *
 *  Created on: 02.05.2012
 *      Author: stephan grein
 */

#include <signal.h>
#include <errno.h>
#include "signal_handler.h"

bool SignalHandler::ExitSignal = false;

SignalHandler::SignalHandler(){ }

SignalHandler::~SignalHandler() { }

const bool SignalHandler::getExitSignal() {
    return ExitSignal;
}

void SignalHandler::setExitSignal() {
    SignalHandler::ExitSignal = ExitSignal;
}

void SignalHandler::setSignalHandler(int ignored) {
    ExitSignal = true;
}

/**
 * Catch and handle CTRL-C
 */
void SignalHandler::setupSignalHandlers() {
    if (signal((int) SIGINT, SignalHandler::setSignalHandler) == SIG_ERR)
        throw SignalException("Error setting up signal handlers.");
}



