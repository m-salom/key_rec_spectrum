#pragma once

poly_t normalize(poly_t h);
void reconstruct_partial(poly_t h, poly_t h_normalized);
void reconstruct_soft(long int nb_dec, double sigma, poly_t h_normalized);


static const double thresholds_200k[12] = {0, 0, 0, -5., -10., -17., -25., -35., -45., -55., -70., -75}; // -N 200000 
static const double thresholds_300k[12] = {0, 0, 0, -3., -10., -12., -20., -25., -40., -42., -55., -60}; // -N 300000 
static const double thresholds_200k_sigma10[12] = {0, 0, 0, -10., -15., -20., -30., -40., -50., -60., -75., -80.}; // -N 200000 -n 10
static const double thresholds_280k_sigma30[12] = {0, 0, 0, -5., -10., -15., -25., -30., -45., -50., -60., -70.}; // -N 280000 -n 30
static const double thresholds_400k_sigma50[12] = {0, 0, 0, -5., -15., -20., -25., -35., -50., -60., -75., -80.}; // -N 400000 -n 50
