#Julian date python module

The juliandate module is designed to simplify the complexities in performing mathematical operations on different dates, while seamlessly interfacing with python's own [datetime](https://docs.python.org/3/library/datetime.html) module.

The juliandate module can be fully functional* on it's own if you wish to only be dependent on a single module. It can also interface with some of the objects found in python's standard datetime module if your project uses them or you wish to take advantage of their functionality as well.

##Aware and Naive
According to the python documentation: 
>Date and time objects may be categorized as "aware" or "naive" depending on whether or not they include timezone information.

Using these definitions the `JulianDate` object can be either aware or naive. If an object is constructed with the timezone set, the object will be aware, otherwise the object is naive.

The one ambiguity here is when a `JulianDate` object is purposely constructed to be in the UTC timezone. While the object is aware of its timezone, internally there is no difference between a timezone of +0 and a naive object. However the user can be certain that any computations done with the former object are accurate and valid.

###Note on time zones:
Unless stated otherwise, a time zone parameter refers to either an integer or floating point number representing the current time zone offset from UTC (including DST).

##Constants
There is a single constant from the juliandate module which represents the J2000 epoch:

juliandate.J2000:  
The J2000 epoch equivalent to `JulianDate(1, 1, 2000, 12, 0, 0, timeZone=0)`

##Types
class juliandate.JulianDate:  
The JulianDate class represents an instance in time and is composed of a day number and fraction of the next day.

##Methods
There is a single module level method which allows the creation of a `JulianDate` object with the current time and 

##`JulianDate` Objects
