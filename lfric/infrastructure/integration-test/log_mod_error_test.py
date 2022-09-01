#!/usr/bin/env python3
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
'''
The Fortran logging module terminates on error. This cannot be tested by the
unit testing framework as it terminates the unit tests as well.
'''
import datetime
import re

from testframework import LFRicLoggingTest, MpiTest, TestEngine, TestFailed

##############################################################################
class log_mod_error_serial_test( MpiTest ):
  '''
  Tests that logging an error terminates execution when run serially.
  '''
  def __init__( self ):
    super(log_mod_error_serial_test, self).__init__( processes=1 )

    self._minimumTimestamp = datetime.datetime.utcnow()

  def test( self, returncode, out, err ):
    expectedLevel = 'ERROR'
    expectedMessage = ' An error was logged.'

    if returncode == 0:
        raise TestFailed('Logging an error did not cause termination to end')
    elif returncode == 127:
        raise TestFailed('Test executable not found')
    elif returncode > 128:
        raise TestFailed('Execution fault such as segmentation fault')

    try:
      timestampString, level, report = err.split( ':', 2 )
      timestampWithoutTimezone = timestampString[:-5]

      timestamp = datetime.datetime.strptime( timestampWithoutTimezone, \
                                              '%Y%m%d%H%M%S.%f' )
    except Exception as ex:
      raise TestFailed( f"Unable to get timestamp from message: {err}", ex )

    if timestamp < self._minimumTimestamp:
      message = 'Expected a timestamp after {} but read {}'
      raise TestFailed( message.format( minimumTimestamp, timestamp ) )

    if level != expectedLevel:
      message = 'Expected "{}" but read "{}"'
      raise TestFailed( message.format( expectedLevel, level ) )

    # We only check the first line as compilers tend to print the return code
    # as well. This will remain true until we can use Fortran 2008 and
    # "stop error".
    #
    first, newline, rest = report.partition( '\n' )
    if first != expectedMessage:
      message = 'Expected "{}" but read "{}"'
      raise TestFailed( message.format( expectedMessage, first ) )

    message = 'Logging an error caused exit as expected with code {code}'
    return message.format( code=returncode )

##############################################################################
class log_mod_error_parallel_test( LFRicLoggingTest ):
  '''
  Tests that logging an error terminates execution when run in parallel.
  '''
  def __init__( self ):
    super(log_mod_error_parallel_test, self).__init__( processes=2 )

    self._minimumTimestamp = datetime.datetime.now(datetime.timezone.utc)
    self._linePattern = re.compile( r'(\d{4})(\d\d)(\d\d)(\d\d)(\d\d)(\d\d)\.(\d{3})([+-])(\d\d)(\d\d):P(\d+):\s*(\w+):\s+(.+)' )

  def test( self, returncode, out, err ):
    expectedLevel = 'ERROR'
    expectedMessage = 'An error was logged.'

    if returncode == 0:
      raise TestFailed( 'Logging an error did not cause termination to end' )

    if out != '':
      message = 'Expected no output on standard out:\n' \
                 + 'Standard out: {out}'
      raise TestFailed( message.format( out=out, err=err ) )

    # We remove the first line as it will be a spin-up message.
    petLog = self.getLFRicLoggingLog()
    petLog = '\n'.join( petLog.splitlines()[1:] )

    match = self._linePattern.match( petLog )
    if match:
      try:
        tzsign = -1 if match.group(8) == '-' else 1
        tzhours = int(match.group(9))
        tzmins = int(match.group(10))
        timezone = datetime.timezone(tzsign
                                     * datetime.timedelta(hours=tzhours,
                                                          minutes=tzmins))
        timestamp = datetime.datetime( int(match.group(1)), # Year
                                       int(match.group(2)), # Month
                                       int(match.group(3)), # Day
                                       int(match.group(4)), # Hour
                                       int(match.group(5)), # Minute
                                       int(match.group(6)), # Second
                                    int(match.group(7)) * 1000, # Microseconds
                                       timezone ) # Timezone
      except Exception as ex:
        raise TestFailed( 'Bad timestamp format: {}'.format( petLog ), ex )
      process = int(match.group(11))
      level = match.group(12)
      report = match.group(13)
    else:
      raise TestFailed( 'Unexpected log message: {}'.format( petLog ) )

    if timestamp < self._minimumTimestamp:
      message = 'Expected a timestamp after {} but read {}'
      raise TestFailed( message.format( minimumTimestamp, timestamp ) )

    if process < 0:
      message = 'Process number went negative'
      raise TestFailed( message )

    if level != expectedLevel:
      message = 'Expected "{}" but read "{}"'
      raise TestFailed( message.format( expectedLevel, level ) )

    # We only check the first line as compilers tend to print the return code
    # as well. This will remain true until we can use Fortran 2008 and
    # "stop error".
    #
    first, newline, rest = report.partition( '\n' )
    if first != expectedMessage:
      message = 'Expected "{}" but read "{}"'
      raise TestFailed( message.format( expectedMessage, first ) )

    message = 'Logging an error caused exit as expected with code {code}'
    return message.format( code=returncode )

##############################################################################
if __name__ == '__main__':
  TestEngine.run( log_mod_error_serial_test() )
  TestEngine.run( log_mod_error_parallel_test() )
