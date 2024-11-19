from django.urls import path
from . import views

urlpatterns = [
    path('', views.home, name='home'),
    path('bisection/', views.bisection, name='bisection'),
    path('regla_falsa/', views.reglaFalsaView, name='regla_falsa'),
    path('PuntoFijo/', views.puntoFijoView, name='PuntoFijo'),
    path('Newton/', views.newtonView, name='Newton'),
    path('secante/', views.secante, name='secante'),
    path('RaicesMultiples/', views.raicesMultiplesView, name='RaicesMultiples'),
    path('iterativos/',views.iterativos,name="iterativos"),
    path('Interpolacion/',views.interpolacion,name="interpolacion"),
    path('lu_methods/', views.lu_methods, name='lu_methods'),
    path('incremental_search/', views.incremental_search,name='incremental_search'),
    path('gaussian_methods/', views.gaussian_elimination_methods, name='gaussian_methods'),

]
