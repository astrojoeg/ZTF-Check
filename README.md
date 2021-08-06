# Welcome to ZTF Check!

With this package you can quickly query the ZTF and Pan-STARRS (PS1) databases for quick reference images, the latter of which is multicolor. All you need is your target's R.A. and declination!

ZTF Check is especially useful when checking to see if you target has an extremely close neighbor or neighboring bright star that could be contaminating it's ZTF light curves.

To run ZTF Check, operate the following terminal command:   
```$ python ztfcheck.py -ra <ra in degrees> -dec <dec in degrees> -q <whether to show Pan-STARRS query results>```

Note that this package makes use of some functions developed by [Zach Vanderbosch](https://github.com/zvanderbosch).

This package is a fun side project and largely a work-in-progress. So, if anything does not work or you have any functionality requests, please correspond with me at jaguidry[at]bu.[dot].edu
