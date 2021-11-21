# mlat-server

# Changes by 'wiedehopf':
 * Ignore servers with bad clock via karma system:
    * Clock sync is checked every 30s, if receiver has bad sync, gets penalty added to its score
    * Badness ranges from "0" to "6.0".  (6 = 30 min timeout)
    * Each time sync is good, bad karma is decremented by small amount.


# Changes by 'TanerH':
1) in `coordinator.py`, write data to a temp file before overwriting the old file (prevents race conditions)
2) in `jsonclient.py`, sanity check the username that is sent, disconnecting if invalid (non-alnum)
3) Clean up exception propagation (dummy loop exception handler added), and clean up the "cleanup" when clients exit.

## Original README below --

This is a Mode S multilateration server that is designed to operate with
clients that do _not_ have synchronized clocks.

It uses ADS-B aircraft that are transmitting DF17 extended squitter position
messages as reference beacons and uses the relative arrival times of those
messages to model the clock characteristics of each receiver.

Then it does multilateration of aircraft that are transmitting only Mode S
using the same receivers.

## License

It is important that you read this section before using or modifying the server!

The server code is licensed under the Affero GPL v3. This license is similar
to the GPL v3, but it has an additional requirement that you must provide
source code to _users who access the server over a network_.

So if you are planning to operate a copy of this server, you must release any
modifications you make to the source code to your users, even if you wouldn't
normally distribute it.

If you are not willing to distribute your changes, you have three options:

 * Contact the copyright holder (Oliver) to discuss a separate license for
   the server code; or
 * Don't allow anyone else to connect to your server, i.e. run only your
   own receivers; or
 * Don't use this server as a basis for your work at all.

The server will automatically provide details of the AGPL license and a link
to the server code, to each client that connects. This is configured in
mlat/config.py. If you make modifications, the suggested process is:

 * Put the modified source code somewhere public (github may be simplest).
 * Update the URL configured in mlat/config.py to point to your modified code.

None of this requires that you make your server publically accessible. If you
want to run a private server with a closed user group, that's fine. But you
must still make the source code for your modified server available to your
users, and they may redistribute it further if they wish.

## Prerequisites

 * Python 3.4 or later. You need the asyncio module which was introduced in 3.4.
 * Numpy and Scipy
 * pygraph (https://github.com/pmatiello/python-graph)
 * pykalman (https://github.com/pykalman/pykalman)
 * optionally, objgraph (https://mg.pov.lt/objgraph/) for leak checking
 * gcc
 * uvloop, ujson, Cython

## Example of how to make it run with virtualenv:

```
apt install python3-pip python3 python3-venv gcc
VENV=/opt/mlat-python-venv
rm -rf $VENV
python3 -m venv $VENV
source $VENV/bin/activate
pip3 install -U pip
pip3 install numpy scipy pykalman python-graph-core uvloop ujson Cython
```

After every code update, recompile the Cython stuff:
```
source $VENV/bin/activate
cd /opt/mlat-server
python3 setup.py build_ext --inplace
```

Starting mlat server:
```
$VENV/bin/python3 /opt/mlat-server/mlat-server
```
(example has git directory cloned into /opt/mlat-server)

For an example service file see systemd-service.example

## Developer-ware

It's all poorly documented and you need to understand quite a bit of the
underlying mathematics of multilateration to make sense of it. Don't expect
to just fire this up and have it all work perfectly first time. You will have
to hack on the code.

## Running

    $ mlat-server --help

## Clients

You need a bunch of receivers running mlat-client:
https://github.com/adsbxchange/mlat-client
The original version by mutability will also work but the adsbexchange client has some changes that are useful.
(https://github.com/mutability/mlat-client)

## Output

Results get passed back to the clients that contributed to the positions.
You can also emit all positions to a local feed, see the command-line help.
