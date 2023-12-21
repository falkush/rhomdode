#define _USE_MATH_DEFINES

#include <SFML/Graphics.hpp>
#include "windows.h" 

void cudathingy(uint8_t* pixels, double pos0, double pos1, double pos2, double vec0, double vec1, double vec2, double addy0, double addy1, double addy2, double addz0, double addz1, double addz2, int remidx, int addidx, int buildidx, uint8_t col, int nbframe);
void cudaInit(bool* blockstmp);
void cudaExit();
void Init3(double* norm0,double* norm1,double* norm2 ,double* point0,double* point1,double* point2,int* mir);

int main()
{
	int i;

	ShowWindow(GetConsoleWindow(), SW_HIDE);
	//ShowWindow(GetConsoleWindow(), SW_SHOW);
	double dist = 1;
	double sqsz = 0.01 / 4;
	double speed = 0.01;
	int nbblocks = 1000;
	int nbframe = 0;


	int mousx, mousy, centralx, centraly,j;

	int remidx;
	double ttmp;
	double distcal;
	int coll;


	double pos0 = 0, pos1 = 500.1, pos2 = 500;
	double vec0, vec1, vec2;
	double addy0, addy1, addy2;
	double addz0, addz1, addz2;
	double multy = (1 - 1280) * sqsz / 2;
	double multz = (720 - 1) * sqsz / 2;

	double buildist = 2;

	double px, py, pz;

	double x[3][3]{};
	double newx[3][3]{};
	x[0][0] = 1;
	x[1][1] = 1;
	x[2][2] = 1;

	for (j = 0; j < 3; j++)
	{
		newx[0][j] = x[0][j] * cos(3*M_PI/2) - sin(3*M_PI / 2) * x[2][j];
		newx[2][j] = x[2][j] * cos(3*M_PI / 2) + sin(3*M_PI / 2) * x[0][j];
		x[0][j] = newx[0][j];
		x[2][j] = newx[2][j];
	}

	
	double anglex, angley;


	bool focus = true;

	sf::RenderWindow window(sf::VideoMode(1280, 720, 32), "Rhombic Dodecahedrons - Press ESC to stop", sf::Style::Titlebar | sf::Style::Close);
	sf::Texture texture;
	sf::Sprite sprite;
	sf::Uint8* pixels = new sf::Uint8[1280 * 720 * 4];
	sf::Vector2i winpos;

	bool* blocks = new bool[nbblocks * nbblocks * nbblocks * sizeof(bool)];
	double* norm0 = new double[12];
	double* norm1 = new double[12];
	double* norm2 = new double[12];
	double* point0 = new double[12];
	double* point1 = new double[12];
	double* point2 = new double[12];
	int* mir = new int[12];

	double qa, qb, qc, discr;
	double t1, t2, tcont;

	int cnx, cny, cnz;

	double cpos0, cpos1, cpos2;
	double tpos0, tpos1, tpos2;
	double tmin;
	int cface;
	bool currblock;
	int blockidx;
	int cnx2;

	bool build = false;
	bool buildfree;
	bool buildset=false;
	bool addblock;

	int buildidx=-1;

	int addidx;

	bool mmm = false;

	int colr=5, colg=0, colb=0;
	uint8_t col;

	cudaInit(blocks);
	Init3(norm0, norm1, norm2, point0, point1, point2, mir);
	
	texture.create(1280, 720);
	window.setMouseCursorVisible(false);

	winpos = window.getPosition();
	SetCursorPos(winpos.x + 1280 / 2, winpos.y + 720 / 2);

	while (window.isOpen())
	{
		//Sleep(1);
		sf::Event event;
		remidx = -1;
		addidx = -1;
		addblock = false;

		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed) window.close();

			if (focus && event.type == sf::Event::MouseButtonPressed && event.mouseButton.button == sf::Mouse::Left)
			{
				if (!build)
				{
					qa = x[0][0] * x[0][0] + x[0][1] * x[0][1] + x[0][2] * x[0][2];
					qb = 2 * (x[0][0] * pos0 + x[0][1] * pos1 + x[0][2] * pos2) - 1000 * (x[0][0] + x[0][1] + x[0][2]);
					qc = pos0 * pos0 + pos1 * pos1 + pos2 * pos2 - 1000 * (pos0 + pos1 + pos2 - 500);

					discr = qb * qb - 4 * qa * qc;

					if (discr > 0)
					{
						t1 = ((-1) * qb - sqrt(discr)) / (2 * qa);
						t2 = ((-1) * qb + sqrt(discr)) / (2 * qa);

						if (!(t1 < 0 && t2 < 0))
						{

							if (t1 * t2 > 0)
							{
								if (t1 < t2) tcont = t1;
								else tcont = t2;

								cpos0 = pos0 + tcont * x[0][0];
								cpos1 = pos1 + tcont * x[0][1];
								cpos2 = pos2 + tcont * x[0][2];
							}
							else
							{
								cpos0 = pos0;
								cpos1 = pos1;
								cpos2 = pos2;
							}

							px = fmod(cpos0, (double)1);
							py = fmod(cpos1, (double)1);
							pz = fmod(cpos2, (double)1);

							if (px < 0) px++;
							if (py < 0) py++;
							if (pz < 0) pz++;

							cnx = cpos0 - px;
							cny = cpos1 - py;
							cnz = cpos2 - pz;


							if ((cnx + cny + cnz) % 2 != 0)
							{
								if (px < py)
								{
									if (px < pz)
									{
										if (px < 1 - py)
										{
											if (px < 1 - pz)
											{
												px++;
												cnx--;
											}
											else
											{
												pz--;
												cnz++;
											}
										}
										else
										{
											if (1 - py < 1 - pz)
											{
												py--;
												cny++;
											}
											else
											{
												pz--;
												cnz++;
											}
										}
									}
									else
									{
										if (pz < 1 - py)
										{
											pz++;
											cnz--;
										}
										else
										{
											py--;
											cny++;
										}
									}
								}
								else
								{
									if (py < pz)
									{
										if (py < 1 - px)
										{
											if (py < 1 - pz)
											{
												py++;
												cny--;
											}
											else
											{
												pz--;
												cnz++;
											}
										}
										else
										{
											if (1 - px < 1 - pz)
											{
												px--;
												cnx++;
											}
											else
											{
												pz--;
												cnz++;
											}
										}
									}
									else
									{
										if (pz < 1 - px)
										{
											pz++;
											cnz--;
										}
										else
										{
											px--;
											cnx++;
										}
									}
								}
							}

							tpos0 = cpos0;
							tpos1 = cpos1;
							tpos2 = cpos2;

							cpos0 = px;
							cpos1 = py;
							cpos2 = pz;

							tmin = 4;


							for (i = 0; i < 12; i++)
							{
								ttmp = (point0[i] - cpos0) * norm0[i] + (point1[i] - cpos1) * norm1[i] + (point2[i] - cpos2) * norm2[i];
								ttmp /= norm0[i] * x[0][0] + norm1[i] * x[0][1] + norm2[i] * x[0][2];

								if (ttmp > 0 && ttmp < tmin)
								{
									tmin = ttmp;
									coll = i;
								}
							}

							cnx -= norm0[coll];
							cny -= norm1[coll];
							cnz -= norm2[coll];

							cpos0 += tmin * x[0][0];
							cpos1 += tmin * x[0][1];
							cpos2 += tmin * x[0][2];

							tpos0 += tmin * x[0][0];
							tpos1 += tmin * x[0][1];
							tpos2 += tmin * x[0][2];

							cpos0 += norm0[coll];
							cpos1 += norm1[coll];
							cpos2 += norm2[coll];

							cface = mir[coll];

							if (cnx < 0 || cnx >= nbblocks || cny < 0 || cny >= nbblocks || cnz < 0 || cnz >= nbblocks)
							{
								currblock = 0;
							}
							else
							{
								if (cnx % 2 != 0) cnx2 = (cnx - 1) / 2;
								else cnx2 = cnx / 2;

								blockidx = cnx2 + nbblocks * cny + nbblocks * nbblocks * cnz;

								currblock = blocks[blockidx];
							}

							while (!currblock && sqrt((tpos0 - 500) * (tpos0 - 500) + (tpos1 - 500) * (tpos1 - 500) + (tpos2 - 500) * (tpos2 - 500)) < 500)
							{
								tmin = 4;


								for (i = 0; i < 12; i++)
								{
									if (i != cface) {
										ttmp = (point0[i] - cpos0) * norm0[i] + (point1[i] - cpos1) * norm1[i] + (point2[i] - cpos2) * norm2[i];
										ttmp /= norm0[i] * x[0][0] + norm1[i] * x[0][1] + norm2[i] * x[0][2];

										if (ttmp > 0 && ttmp < tmin)
										{
											tmin = ttmp;
											coll = i;
										}
									}
								}

								cnx -= norm0[coll];
								cny -= norm1[coll];
								cnz -= norm2[coll];

								cpos0 += tmin * x[0][0];
								cpos1 += tmin * x[0][1];
								cpos2 += tmin * x[0][2];

								tpos0 += tmin * x[0][0];
								tpos1 += tmin * x[0][1];
								tpos2 += tmin * x[0][2];

								cpos0 += norm0[coll];
								cpos1 += norm1[coll];
								cpos2 += norm2[coll];

								cface = mir[coll];

								if (cnx < 0 || cnx >= nbblocks || cny < 0 || cny >= nbblocks || cnz < 0 || cnz >= nbblocks)
								{
									currblock = false;
								}
							    else
								{
									if (cnx % 2 != 0) cnx2 = (cnx - 1) / 2;
									else cnx2 = cnx / 2;
								cnx2 = cnx / 2;

								blockidx = cnx2 + nbblocks * cny + nbblocks * nbblocks * cnz;

								currblock = blocks[blockidx];
								}
							}

							if (currblock)
							{
								remidx = blockidx;
								blocks[remidx] = false;
							}
						}
					}
				}
				else
				{
					if (buildset)
					{
						buildset = false;
						addidx = buildidx;
						buildidx = -1;
						blocks[addidx] = true;
						addblock = true;
					}
				}
			}
			if (focus && event.type == sf::Event::MouseMoved)
			{
				POINT p;
				GetCursorPos(&p);
				winpos = window.getPosition();
				centralx = winpos.x + 1280 / 2;
				centraly = winpos.y + 720 / 2;
				SetCursorPos(centralx, centraly);

				mousx = p.x - centralx;
				mousy = p.y - centraly;

				anglex = 0.002 * mousx;
				angley = 0.002 * mousy;

				if (mousx > 0)
				{
					for (j = 0; j < 3; j++)
					{
						newx[0][j] = x[0][j] * cos(anglex) + sin(anglex) * x[1][j];
						newx[1][j] = x[1][j] * cos(anglex) - sin(anglex) * x[0][j];
						x[0][j] = newx[0][j];
						x[1][j] = newx[1][j];
					}
				}
				else if (mousx < 0)
				{
					anglex = -anglex;

					for (j = 0; j < 3; j++)
					{
						newx[0][j] = x[0][j] * cos(anglex) - sin(anglex) * x[1][j];
						newx[1][j] = x[1][j] * cos(anglex) + sin(anglex) * x[0][j];
						x[0][j] = newx[0][j];
						x[1][j] = newx[1][j];
					}
					
				}

				if (mousy < 0)
				{
					angley = -angley;
					for (j = 0; j < 3; j++)
					{
						newx[0][j] = x[0][j] * cos(angley) + sin(angley) * x[2][j];
						newx[2][j] = x[2][j] * cos(angley) - sin(angley) * x[0][j];
						x[0][j] = newx[0][j];
						x[2][j] = newx[2][j];
					}
					
				}
				else if (mousy > 0)
				{
					for (j = 0; j < 3; j++)
					{
						newx[0][j] = x[0][j] * cos(angley) - sin(angley) * x[2][j];
						newx[2][j] = x[2][j] * cos(angley) + sin(angley) * x[0][j];
						x[0][j] = newx[0][j];
						x[2][j] = newx[2][j];
					}
				}


			}
			if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::M)
			{
				if (speed == 1)speed = 0.01;
				else speed = 1;
			}
			if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::N)
			{
				build = !build;
			}
			if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::Z)
			{
				colr++;
				if (mmm)
				{
					if (colr == 7) colr = 0;
				}
				else
				{
					if (colr == 6) colr = 0;
				}

			}
			if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::X)
			{
				colg++;
				if (colg == 6) colg = 0;
			}
			if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::C)
			{
				colb++;
				if (colb == 6) colb = 0;
			}
		}

		if (sf::Keyboard::isKeyPressed(sf::Keyboard::P))
		{
			
			pos0 = 0; pos1 = 500.1; pos2 = 500;

			x[0][0] = 1; x[0][1] = 0; x[0][2] = 0;
			x[1][0] = 0; x[1][1] = 1; x[1][2] = 0;
			x[2][0] = 0; x[2][1] = 0; x[2][2] = 1;

			for (j = 0; j < 3; j++)
			{
				newx[0][j] = x[0][j] * cos(3 * M_PI / 2) - sin(3 * M_PI / 2) * x[2][j];
				newx[2][j] = x[2][j] * cos(3 * M_PI / 2) + sin(3 * M_PI / 2) * x[0][j];
				x[0][j] = newx[0][j];
				x[2][j] = newx[2][j];
			}

			speed = 0.01;
		}

		if (sf::Keyboard::isKeyPressed(sf::Keyboard::Escape))
		{
			focus = false;
			window.setMouseCursorVisible(true);
		}

		if (sf::Mouse::isButtonPressed(sf::Mouse::Left))
		{
				focus = true;
				window.setMouseCursorVisible(false);
		}



		if (sf::Keyboard::isKeyPressed(sf::Keyboard::W))
		{

			distcal = 0;
			qa = x[0][0] * x[0][0] + x[0][1] * x[0][1] + x[0][2] * x[0][2];
			qb = 2 * (x[0][0] * pos0 + x[0][1] * pos1 + x[0][2] * pos2) - 1000 * (x[0][0] + x[0][1] + x[0][2]);
			qc = pos0 * pos0 + pos1 * pos1 + pos2 * pos2 - 1000 * (pos0 + pos1 + pos2 - 500);

			discr = qb * qb - 4 * qa * qc;

			if (discr <= 0)
			{
				pos0 += x[0][0] * speed;
				pos1 += x[0][1] * speed;
				pos2 += x[0][2] * speed;
			}
			else
			{
				t1 = ((-1) * qb - sqrt(discr)) / (2 * qa);
				t2 = ((-1) * qb + sqrt(discr)) / (2 * qa);

				if (t1 < 0 && t2 < 0)
				{
					pos0 += x[0][0] * speed;
					pos1 += x[0][1] * speed;
					pos2 += x[0][2] * speed;
				}
				else
				{
					if (t1 * t2 > 0)
					{
						if (t1 < t2) tcont = t1;
						else tcont = t2;

						cpos0 = pos0 + tcont * x[0][0];
						cpos1 = pos1 + tcont * x[0][1];
						cpos2 = pos2 + tcont * x[0][2];
						distcal = tcont;
					}
					else
					{
						cpos0 = pos0;
						cpos1 = pos1;
						cpos2 = pos2;
					}

					px = fmod(cpos0, (double)1);
					py = fmod(cpos1, (double)1);
					pz = fmod(cpos2, (double)1);

					if (px < 0) px++;
					if (py < 0) py++;
					if (pz < 0) pz++;

					cnx = cpos0 - px;
					cny = cpos1 - py;
					cnz = cpos2 - pz;


					if ((cnx + cny + cnz) % 2 != 0)
					{
						if (px < py)
						{
							if (px < pz)
							{
								if (px < 1 - py)
								{
									if (px < 1 - pz)
									{
										px++;
										cnx--;
									}
									else
									{
										pz--;
										cnz++;
									}
								}
								else
								{
									if (1 - py < 1 - pz)
									{
										py--;
										cny++;
									}
									else
									{
										pz--;
										cnz++;
									}
								}
							}
							else
							{
								if (pz < 1 - py)
								{
									pz++;
									cnz--;
								}
								else
								{
									py--;
									cny++;
								}
							}
						}
						else
						{
							if (py < pz)
							{
								if (py < 1 - px)
								{
									if (py < 1 - pz)
									{
										py++;
										cny--;
									}
									else
									{
										pz--;
										cnz++;
									}
								}
								else
								{
									if (1 - px < 1 - pz)
									{
										px--;
										cnx++;
									}
									else
									{
										pz--;
										cnz++;
									}
								}
							}
							else
							{
								if (pz < 1 - px)
								{
									pz++;
									cnz--;
								}
								else
								{
									px--;
									cnx++;
								}
							}
						}
					}

					tpos0 = cpos0;
					tpos1 = cpos1;
					tpos2 = cpos2;

					cpos0 = px;
					cpos1 = py;
					cpos2 = pz;

					tmin = 4;


					for (i = 0; i < 12; i++)
					{
						ttmp = (point0[i] - cpos0) * norm0[i] + (point1[i] - cpos1) * norm1[i] + (point2[i] - cpos2) * norm2[i];
						ttmp /= norm0[i] * x[0][0] + norm1[i] * x[0][1] + norm2[i] * x[0][2];

						if (ttmp > 0 && ttmp < tmin)
						{
							tmin = ttmp;
							coll = i;
						}
					}

					cnx -= norm0[coll];
					cny -= norm1[coll];
					cnz -= norm2[coll];

					cpos0 += tmin * x[0][0];
					cpos1 += tmin * x[0][1];
					cpos2 += tmin * x[0][2];

					tpos0 += tmin * x[0][0];
					tpos1 += tmin * x[0][1];
					tpos2 += tmin * x[0][2];

					cpos0 += norm0[coll];
					cpos1 += norm1[coll];
					cpos2 += norm2[coll];

					distcal += tmin;

					cface = mir[coll];

					if (cnx < 0 || cnx >= nbblocks || cny < 0 || cny >= nbblocks || cnz < 0 || cnz >= nbblocks)
					{
						currblock = 0;
					}
					else
					{
						if (cnx % 2 != 0) cnx2 = (cnx - 1) / 2;
						else cnx2 = cnx / 2;

						blockidx = cnx2 + nbblocks * cny + nbblocks * nbblocks * cnz;

						currblock = blocks[blockidx];
					}

					while (!currblock && sqrt((tpos0 - 500) * (tpos0 - 500) + (tpos1 - 500) * (tpos1 - 500) + (tpos2 - 500) * (tpos2 - 500)) < 500)
					{
						tmin = 4;


						for (i = 0; i < 12; i++)
						{
							if (i != cface) {
								ttmp = (point0[i] - cpos0) * norm0[i] + (point1[i] - cpos1) * norm1[i] + (point2[i] - cpos2) * norm2[i];
								ttmp /= norm0[i] * x[0][0] + norm1[i] * x[0][1] + norm2[i] * x[0][2];

								if (ttmp > 0 && ttmp < tmin)
								{
									tmin = ttmp;
									coll = i;
								}
							}
						}

						cnx -= norm0[coll];
						cny -= norm1[coll];
						cnz -= norm2[coll];

						cpos0 += tmin * x[0][0];
						cpos1 += tmin * x[0][1];
						cpos2 += tmin * x[0][2];

						tpos0 += tmin * x[0][0];
						tpos1 += tmin * x[0][1];
						tpos2 += tmin * x[0][2];

						distcal += tmin;

						cpos0 += norm0[coll];
						cpos1 += norm1[coll];
						cpos2 += norm2[coll];

						cface = mir[coll];

						if (cnx < 0 || cnx >= nbblocks || cny < 0 || cny >= nbblocks || cnz < 0 || cnz >= nbblocks)
						{
							currblock = false;
						}
						else
						{
							if (cnx % 2 != 0) cnx2 = (cnx - 1) / 2;
							else cnx2 = cnx / 2;
						cnx2 = cnx / 2;

						blockidx = cnx2 + nbblocks * cny + nbblocks * nbblocks * cnz;

						currblock = blocks[blockidx];
						}
					}

					if (currblock)
					{
						if (speed < distcal)
						{
							pos0 += x[0][0] * speed;
							pos1 += x[0][1] * speed;
							pos2 += x[0][2] * speed;
						}
					}
					else
					{
						pos0 += x[0][0] * speed;
						pos1 += x[0][1] * speed;
						pos2 += x[0][2] * speed;
					}
				}
			}
		}

		if (sf::Keyboard::isKeyPressed(sf::Keyboard::D))
		{
			newx[0][0] = x[0][0];
			newx[0][1] = x[0][1];
			newx[0][2] = x[0][2];

			x[0][0] = x[1][0];
			x[0][1] = x[1][1];
			x[0][2] = x[1][2];

			distcal = 0;
			qa = x[0][0] * x[0][0] + x[0][1] * x[0][1] + x[0][2] * x[0][2];
			qb = 2 * (x[0][0] * pos0 + x[0][1] * pos1 + x[0][2] * pos2) - 1000 * (x[0][0] + x[0][1] + x[0][2]);
			qc = pos0 * pos0 + pos1 * pos1 + pos2 * pos2 - 1000 * (pos0 + pos1 + pos2 - 500);

			discr = qb * qb - 4 * qa * qc;

			if (discr <= 0)
			{
				pos0 += x[0][0] * speed;
				pos1 += x[0][1] * speed;
				pos2 += x[0][2] * speed;
			}
			else
			{
				t1 = ((-1) * qb - sqrt(discr)) / (2 * qa);
				t2 = ((-1) * qb + sqrt(discr)) / (2 * qa);

				if (t1 < 0 && t2 < 0)
				{
					pos0 += x[0][0] * speed;
					pos1 += x[0][1] * speed;
					pos2 += x[0][2] * speed;
				}
				else
				{
					if (t1 * t2 > 0)
					{
						if (t1 < t2) tcont = t1;
						else tcont = t2;

						cpos0 = pos0 + tcont * x[0][0];
						cpos1 = pos1 + tcont * x[0][1];
						cpos2 = pos2 + tcont * x[0][2];
						distcal = tcont;
					}
					else
					{
						cpos0 = pos0;
						cpos1 = pos1;
						cpos2 = pos2;
					}

					px = fmod(cpos0, (double)1);
					py = fmod(cpos1, (double)1);
					pz = fmod(cpos2, (double)1);

					if (px < 0) px++;
					if (py < 0) py++;
					if (pz < 0) pz++;

					cnx = cpos0 - px;
					cny = cpos1 - py;
					cnz = cpos2 - pz;


					if ((cnx + cny + cnz) % 2 != 0)
					{
						if (px < py)
						{
							if (px < pz)
							{
								if (px < 1 - py)
								{
									if (px < 1 - pz)
									{
										px++;
										cnx--;
									}
									else
									{
										pz--;
										cnz++;
									}
								}
								else
								{
									if (1 - py < 1 - pz)
									{
										py--;
										cny++;
									}
									else
									{
										pz--;
										cnz++;
									}
								}
							}
							else
							{
								if (pz < 1 - py)
								{
									pz++;
									cnz--;
								}
								else
								{
									py--;
									cny++;
								}
							}
						}
						else
						{
							if (py < pz)
							{
								if (py < 1 - px)
								{
									if (py < 1 - pz)
									{
										py++;
										cny--;
									}
									else
									{
										pz--;
										cnz++;
									}
								}
								else
								{
									if (1 - px < 1 - pz)
									{
										px--;
										cnx++;
									}
									else
									{
										pz--;
										cnz++;
									}
								}
							}
							else
							{
								if (pz < 1 - px)
								{
									pz++;
									cnz--;
								}
								else
								{
									px--;
									cnx++;
								}
							}
						}
					}

					tpos0 = cpos0;
					tpos1 = cpos1;
					tpos2 = cpos2;

					cpos0 = px;
					cpos1 = py;
					cpos2 = pz;

					tmin = 4;


					for (i = 0; i < 12; i++)
					{
						ttmp = (point0[i] - cpos0) * norm0[i] + (point1[i] - cpos1) * norm1[i] + (point2[i] - cpos2) * norm2[i];
						ttmp /= norm0[i] * x[0][0] + norm1[i] * x[0][1] + norm2[i] * x[0][2];

						if (ttmp > 0 && ttmp < tmin)
						{
							tmin = ttmp;
							coll = i;
						}
					}

					cnx -= norm0[coll];
					cny -= norm1[coll];
					cnz -= norm2[coll];

					cpos0 += tmin * x[0][0];
					cpos1 += tmin * x[0][1];
					cpos2 += tmin * x[0][2];

					tpos0 += tmin * x[0][0];
					tpos1 += tmin * x[0][1];
					tpos2 += tmin * x[0][2];

					cpos0 += norm0[coll];
					cpos1 += norm1[coll];
					cpos2 += norm2[coll];

					distcal += tmin;

					cface = mir[coll];

					if (cnx < 0 || cnx >= nbblocks || cny < 0 || cny >= nbblocks || cnz < 0 || cnz >= nbblocks)
					{
						currblock = 0;
					}
					else
					{
						if (cnx % 2 != 0) cnx2 = (cnx - 1) / 2;
						else cnx2 = cnx / 2;

						blockidx = cnx2 + nbblocks * cny + nbblocks * nbblocks * cnz;

						currblock = blocks[blockidx];
					}

					while (!currblock && sqrt((tpos0 - 500) * (tpos0 - 500) + (tpos1 - 500) * (tpos1 - 500) + (tpos2 - 500) * (tpos2 - 500)) < 500)
					{
						tmin = 4;


						for (i = 0; i < 12; i++)
						{
							if (i != cface) {
								ttmp = (point0[i] - cpos0) * norm0[i] + (point1[i] - cpos1) * norm1[i] + (point2[i] - cpos2) * norm2[i];
								ttmp /= norm0[i] * x[0][0] + norm1[i] * x[0][1] + norm2[i] * x[0][2];

								if (ttmp > 0 && ttmp < tmin)
								{
									tmin = ttmp;
									coll = i;
								}
							}
						}

						cnx -= norm0[coll];
						cny -= norm1[coll];
						cnz -= norm2[coll];

						cpos0 += tmin * x[0][0];
						cpos1 += tmin * x[0][1];
						cpos2 += tmin * x[0][2];

						tpos0 += tmin * x[0][0];
						tpos1 += tmin * x[0][1];
						tpos2 += tmin * x[0][2];

						distcal += tmin;

						cpos0 += norm0[coll];
						cpos1 += norm1[coll];
						cpos2 += norm2[coll];

						cface = mir[coll];

						if (cnx < 0 || cnx >= nbblocks || cny < 0 || cny >= nbblocks || cnz < 0 || cnz >= nbblocks)
						{
							currblock = false;
						}
						else
						{
							if (cnx % 2 != 0) cnx2 = (cnx - 1) / 2;
							else cnx2 = cnx / 2;
						cnx2 = cnx / 2;

						blockidx = cnx2 + nbblocks * cny + nbblocks * nbblocks * cnz;

						currblock = blocks[blockidx];
						}
					}

					if (currblock)
					{
						if (speed < distcal)
						{
							pos0 += x[0][0] * speed;
							pos1 += x[0][1] * speed;
							pos2 += x[0][2] * speed;
						}
					}
					else
					{
						pos0 += x[0][0] * speed;
						pos1 += x[0][1] * speed;
						pos2 += x[0][2] * speed;
					}
				}
			}

			x[0][0] = newx[0][0];
			x[0][1] = newx[0][1];
			x[0][2] = newx[0][2];

		}

		if (sf::Keyboard::isKeyPressed(sf::Keyboard::A))
		{
			newx[0][0] = x[0][0];
			newx[0][1] = x[0][1];
			newx[0][2] = x[0][2];

			x[0][0] = -x[1][0];
			x[0][1] = -x[1][1];
			x[0][2] = -x[1][2];

			distcal = 0;
			qa = x[0][0] * x[0][0] + x[0][1] * x[0][1] + x[0][2] * x[0][2];
			qb = 2 * (x[0][0] * pos0 + x[0][1] * pos1 + x[0][2] * pos2) - 1000 * (x[0][0] + x[0][1] + x[0][2]);
			qc = pos0 * pos0 + pos1 * pos1 + pos2 * pos2 - 1000 * (pos0 + pos1 + pos2 - 500);

			discr = qb * qb - 4 * qa * qc;

			if (discr <= 0)
			{
				pos0 += x[0][0] * speed;
				pos1 += x[0][1] * speed;
				pos2 += x[0][2] * speed;
			}
			else
			{
				t1 = ((-1) * qb - sqrt(discr)) / (2 * qa);
				t2 = ((-1) * qb + sqrt(discr)) / (2 * qa);

				if (t1 < 0 && t2 < 0)
				{
					pos0 += x[0][0] * speed;
					pos1 += x[0][1] * speed;
					pos2 += x[0][2] * speed;
				}
				else
				{
					if (t1 * t2 > 0)
					{
						if (t1 < t2) tcont = t1;
						else tcont = t2;

						cpos0 = pos0 + tcont * x[0][0];
						cpos1 = pos1 + tcont * x[0][1];
						cpos2 = pos2 + tcont * x[0][2];
						distcal = tcont;
					}
					else
					{
						cpos0 = pos0;
						cpos1 = pos1;
						cpos2 = pos2;
					}

					px = fmod(cpos0, (double)1);
					py = fmod(cpos1, (double)1);
					pz = fmod(cpos2, (double)1);

					if (px < 0) px++;
					if (py < 0) py++;
					if (pz < 0) pz++;

					cnx = cpos0 - px;
					cny = cpos1 - py;
					cnz = cpos2 - pz;


					if ((cnx + cny + cnz) % 2 != 0)
					{
						if (px < py)
						{
							if (px < pz)
							{
								if (px < 1 - py)
								{
									if (px < 1 - pz)
									{
										px++;
										cnx--;
									}
									else
									{
										pz--;
										cnz++;
									}
								}
								else
								{
									if (1 - py < 1 - pz)
									{
										py--;
										cny++;
									}
									else
									{
										pz--;
										cnz++;
									}
								}
							}
							else
							{
								if (pz < 1 - py)
								{
									pz++;
									cnz--;
								}
								else
								{
									py--;
									cny++;
								}
							}
						}
						else
						{
							if (py < pz)
							{
								if (py < 1 - px)
								{
									if (py < 1 - pz)
									{
										py++;
										cny--;
									}
									else
									{
										pz--;
										cnz++;
									}
								}
								else
								{
									if (1 - px < 1 - pz)
									{
										px--;
										cnx++;
									}
									else
									{
										pz--;
										cnz++;
									}
								}
							}
							else
							{
								if (pz < 1 - px)
								{
									pz++;
									cnz--;
								}
								else
								{
									px--;
									cnx++;
								}
							}
						}
					}

					tpos0 = cpos0;
					tpos1 = cpos1;
					tpos2 = cpos2;

					cpos0 = px;
					cpos1 = py;
					cpos2 = pz;

					tmin = 4;


					for (i = 0; i < 12; i++)
					{
						ttmp = (point0[i] - cpos0) * norm0[i] + (point1[i] - cpos1) * norm1[i] + (point2[i] - cpos2) * norm2[i];
						ttmp /= norm0[i] * x[0][0] + norm1[i] * x[0][1] + norm2[i] * x[0][2];

						if (ttmp > 0 && ttmp < tmin)
						{
							tmin = ttmp;
							coll = i;
						}
					}

					cnx -= norm0[coll];
					cny -= norm1[coll];
					cnz -= norm2[coll];

					cpos0 += tmin * x[0][0];
					cpos1 += tmin * x[0][1];
					cpos2 += tmin * x[0][2];

					tpos0 += tmin * x[0][0];
					tpos1 += tmin * x[0][1];
					tpos2 += tmin * x[0][2];

					cpos0 += norm0[coll];
					cpos1 += norm1[coll];
					cpos2 += norm2[coll];

					distcal += tmin;

					cface = mir[coll];

					if (cnx < 0 || cnx >= nbblocks || cny < 0 || cny >= nbblocks || cnz < 0 || cnz >= nbblocks)
					{
						currblock = 0;
					}
					else
					{
						if (cnx % 2 != 0) cnx2 = (cnx - 1) / 2;
						else cnx2 = cnx / 2;

						blockidx = cnx2 + nbblocks * cny + nbblocks * nbblocks * cnz;

						currblock = blocks[blockidx];
					}

					while (!currblock && sqrt((tpos0 - 500) * (tpos0 - 500) + (tpos1 - 500) * (tpos1 - 500) + (tpos2 - 500) * (tpos2 - 500)) < 500)
					{
						tmin = 4;


						for (i = 0; i < 12; i++)
						{
							if (i != cface) {
								ttmp = (point0[i] - cpos0) * norm0[i] + (point1[i] - cpos1) * norm1[i] + (point2[i] - cpos2) * norm2[i];
								ttmp /= norm0[i] * x[0][0] + norm1[i] * x[0][1] + norm2[i] * x[0][2];

								if (ttmp > 0 && ttmp < tmin)
								{
									tmin = ttmp;
									coll = i;
								}
							}
						}

						cnx -= norm0[coll];
						cny -= norm1[coll];
						cnz -= norm2[coll];

						cpos0 += tmin * x[0][0];
						cpos1 += tmin * x[0][1];
						cpos2 += tmin * x[0][2];

						tpos0 += tmin * x[0][0];
						tpos1 += tmin * x[0][1];
						tpos2 += tmin * x[0][2];

						distcal += tmin;

						cpos0 += norm0[coll];
						cpos1 += norm1[coll];
						cpos2 += norm2[coll];

						cface = mir[coll];

						if (cnx < 0 || cnx >= nbblocks || cny < 0 || cny >= nbblocks || cnz < 0 || cnz >= nbblocks)
						{
							currblock = false;
						}
						else
						{
							if (cnx % 2 != 0) cnx2 = (cnx - 1) / 2;
							else cnx2 = cnx / 2;
						cnx2 = cnx / 2;

						blockidx = cnx2 + nbblocks * cny + nbblocks * nbblocks * cnz;

						currblock = blocks[blockidx];
						}
					}

					if (currblock)
					{
						if (speed < distcal)
						{
							pos0 += x[0][0] * speed;
							pos1 += x[0][1] * speed;
							pos2 += x[0][2] * speed;
						}
					}
					else
					{
						pos0 += x[0][0] * speed;
						pos1 += x[0][1] * speed;
						pos2 += x[0][2] * speed;
					}
				}
			}

			x[0][0] = newx[0][0];
			x[0][1] = newx[0][1];
			x[0][2] = newx[0][2];
		}

		if (sf::Keyboard::isKeyPressed(sf::Keyboard::S))
		{
			x[0][0] *= -1;
			x[0][1] *= -1;
			x[0][2] *= -1;

			distcal = 0;
			qa = x[0][0] * x[0][0] + x[0][1] * x[0][1] + x[0][2] * x[0][2];
			qb = 2 * (x[0][0] * pos0 + x[0][1] * pos1 + x[0][2] * pos2) - 1000 * (x[0][0] + x[0][1] + x[0][2]);
			qc = pos0 * pos0 + pos1 * pos1 + pos2 * pos2 - 1000 * (pos0 + pos1 + pos2 - 500);

			discr = qb * qb - 4 * qa * qc;

			if (discr <= 0)
			{
				pos0 += x[0][0] * speed;
				pos1 += x[0][1] * speed;
				pos2 += x[0][2] * speed;
			}
			else
			{
				t1 = ((-1) * qb - sqrt(discr)) / (2 * qa);
				t2 = ((-1) * qb + sqrt(discr)) / (2 * qa);

				if (t1 < 0 && t2 < 0)
				{
					pos0 += x[0][0] * speed;
					pos1 += x[0][1] * speed;
					pos2 += x[0][2] * speed;
				}
				else
				{
					if (t1 * t2 > 0)
					{
						if (t1 < t2) tcont = t1;
						else tcont = t2;

						cpos0 = pos0 + tcont * x[0][0];
						cpos1 = pos1 + tcont * x[0][1];
						cpos2 = pos2 + tcont * x[0][2];
						distcal = tcont;
					}
					else
					{
						cpos0 = pos0;
						cpos1 = pos1;
						cpos2 = pos2;
					}

					px = fmod(cpos0, (double)1);
					py = fmod(cpos1, (double)1);
					pz = fmod(cpos2, (double)1);

					if (px < 0) px++;
					if (py < 0) py++;
					if (pz < 0) pz++;

					cnx = cpos0 - px;
					cny = cpos1 - py;
					cnz = cpos2 - pz;


					if ((cnx + cny + cnz) % 2 != 0)
					{
						if (px < py)
						{
							if (px < pz)
							{
								if (px < 1 - py)
								{
									if (px < 1 - pz)
									{
										px++;
										cnx--;
									}
									else
									{
										pz--;
										cnz++;
									}
								}
								else
								{
									if (1 - py < 1 - pz)
									{
										py--;
										cny++;
									}
									else
									{
										pz--;
										cnz++;
									}
								}
							}
							else
							{
								if (pz < 1 - py)
								{
									pz++;
									cnz--;
								}
								else
								{
									py--;
									cny++;
								}
							}
						}
						else
						{
							if (py < pz)
							{
								if (py < 1 - px)
								{
									if (py < 1 - pz)
									{
										py++;
										cny--;
									}
									else
									{
										pz--;
										cnz++;
									}
								}
								else
								{
									if (1 - px < 1 - pz)
									{
										px--;
										cnx++;
									}
									else
									{
										pz--;
										cnz++;
									}
								}
							}
							else
							{
								if (pz < 1 - px)
								{
									pz++;
									cnz--;
								}
								else
								{
									px--;
									cnx++;
								}
							}
						}
					}

					tpos0 = cpos0;
					tpos1 = cpos1;
					tpos2 = cpos2;

					cpos0 = px;
					cpos1 = py;
					cpos2 = pz;

					tmin = 4;


					for (i = 0; i < 12; i++)
					{
						ttmp = (point0[i] - cpos0) * norm0[i] + (point1[i] - cpos1) * norm1[i] + (point2[i] - cpos2) * norm2[i];
						ttmp /= norm0[i] * x[0][0] + norm1[i] * x[0][1] + norm2[i] * x[0][2];

						if (ttmp > 0 && ttmp < tmin)
						{
							tmin = ttmp;
							coll = i;
						}
					}

					cnx -= norm0[coll];
					cny -= norm1[coll];
					cnz -= norm2[coll];

					cpos0 += tmin * x[0][0];
					cpos1 += tmin * x[0][1];
					cpos2 += tmin * x[0][2];

					tpos0 += tmin * x[0][0];
					tpos1 += tmin * x[0][1];
					tpos2 += tmin * x[0][2];

					cpos0 += norm0[coll];
					cpos1 += norm1[coll];
					cpos2 += norm2[coll];

					distcal += tmin;

					cface = mir[coll];

					if (cnx < 0 || cnx >= nbblocks || cny < 0 || cny >= nbblocks || cnz < 0 || cnz >= nbblocks)
					{
						currblock = 0;
					}
					else
					{
						if (cnx % 2 != 0) cnx2 = (cnx - 1) / 2;
						else cnx2 = cnx / 2;

						blockidx = cnx2 + nbblocks * cny + nbblocks * nbblocks * cnz;

						currblock = blocks[blockidx];
					}

					while (!currblock && sqrt((tpos0 - 500) * (tpos0 - 500) + (tpos1 - 500) * (tpos1 - 500) + (tpos2 - 500) * (tpos2 - 500)) < 500)
					{
						tmin = 4;


						for (i = 0; i < 12; i++)
						{
							if (i != cface) {
								ttmp = (point0[i] - cpos0) * norm0[i] + (point1[i] - cpos1) * norm1[i] + (point2[i] - cpos2) * norm2[i];
								ttmp /= norm0[i] * x[0][0] + norm1[i] * x[0][1] + norm2[i] * x[0][2];

								if (ttmp > 0 && ttmp < tmin)
								{
									tmin = ttmp;
									coll = i;
								}
							}
						}

						cnx -= norm0[coll];
						cny -= norm1[coll];
						cnz -= norm2[coll];

						cpos0 += tmin * x[0][0];
						cpos1 += tmin * x[0][1];
						cpos2 += tmin * x[0][2];

						tpos0 += tmin * x[0][0];
						tpos1 += tmin * x[0][1];
						tpos2 += tmin * x[0][2];

						distcal += tmin;

						cpos0 += norm0[coll];
						cpos1 += norm1[coll];
						cpos2 += norm2[coll];

						cface = mir[coll];

						if (cnx < 0 || cnx >= nbblocks || cny < 0 || cny >= nbblocks || cnz < 0 || cnz >= nbblocks)
						{
							currblock = false;
						}
						else
						{
							if (cnx % 2 != 0) cnx2 = (cnx - 1) / 2;
							else cnx2 = cnx / 2;
						cnx2 = cnx / 2;

						blockidx = cnx2 + nbblocks * cny + nbblocks * nbblocks * cnz;

						currblock = blocks[blockidx];
						}
					}

					if (currblock)
					{
						if (speed < distcal)
						{
							pos0 += x[0][0] * speed;
							pos1 += x[0][1] * speed;
							pos2 += x[0][2] * speed;
						}
					}
					else
					{
						pos0 += x[0][0] * speed;
						pos1 += x[0][1] * speed;
						pos2 += x[0][2] * speed;
					}
				}
			}

			x[0][0] *= -1;
			x[0][1] *= -1;
			x[0][2] *= -1;
		}

		if (sf::Keyboard::isKeyPressed(sf::Keyboard::V))
		{
			if(buildist>2) buildist -= 0.1;
		}
		if (sf::Keyboard::isKeyPressed(sf::Keyboard::B))
		{
			buildist += 0.1;
		}


		if (focus)
		{
			nbframe++;

			vec0 = dist * x[0][0] + multy * x[1][0] + multz * x[2][0];
			vec1 = dist * x[0][1] + multy * x[1][1] + multz * x[2][1];
			vec2 = dist * x[0][2] + multy * x[1][2] + multz * x[2][2];

			addy0 = sqsz * x[1][0];
			addy1 = sqsz * x[1][1];
			addy2 = sqsz * x[1][2];

			addz0 = -sqsz * x[2][0];
			addz1 = -sqsz * x[2][1];
			addz2 = -sqsz * x[2][2];

			if (build)
			{
				if (!addblock)
				{
					buildfree = false;

					distcal = 0;
					qa = x[0][0] * x[0][0] + x[0][1] * x[0][1] + x[0][2] * x[0][2];
					qb = 2 * (x[0][0] * pos0 + x[0][1] * pos1 + x[0][2] * pos2) - 1000 * (x[0][0] + x[0][1] + x[0][2]);
					qc = pos0 * pos0 + pos1 * pos1 + pos2 * pos2 - 1000 * (pos0 + pos1 + pos2 - 500);

					discr = qb * qb - 4 * qa * qc;

					if (discr <= 0)
					{
						buildfree = true;
					}
					else
					{
						t1 = ((-1) * qb - sqrt(discr)) / (2 * qa);
						t2 = ((-1) * qb + sqrt(discr)) / (2 * qa);

						if (t1 < 0 && t2 < 0)
						{
							buildfree = true;
						}
						else
						{
							if (t1 * t2 > 0)
							{
								if (t1 < t2) tcont = t1;
								else tcont = t2;

								cpos0 = pos0 + tcont * x[0][0];
								cpos1 = pos1 + tcont * x[0][1];
								cpos2 = pos2 + tcont * x[0][2];
								distcal = tcont;
							}
							else
							{
								cpos0 = pos0;
								cpos1 = pos1;
								cpos2 = pos2;
							}

							px = fmod(cpos0, (double)1);
							py = fmod(cpos1, (double)1);
							pz = fmod(cpos2, (double)1);

							if (px < 0) px++;
							if (py < 0) py++;
							if (pz < 0) pz++;

							cnx = cpos0 - px;
							cny = cpos1 - py;
							cnz = cpos2 - pz;


							if ((cnx + cny + cnz) % 2 != 0)
							{
								if (px < py)
								{
									if (px < pz)
									{
										if (px < 1 - py)
										{
											if (px < 1 - pz)
											{
												px++;
												cnx--;
											}
											else
											{
												pz--;
												cnz++;
											}
										}
										else
										{
											if (1 - py < 1 - pz)
											{
												py--;
												cny++;
											}
											else
											{
												pz--;
												cnz++;
											}
										}
									}
									else
									{
										if (pz < 1 - py)
										{
											pz++;
											cnz--;
										}
										else
										{
											py--;
											cny++;
										}
									}
								}
								else
								{
									if (py < pz)
									{
										if (py < 1 - px)
										{
											if (py < 1 - pz)
											{
												py++;
												cny--;
											}
											else
											{
												pz--;
												cnz++;
											}
										}
										else
										{
											if (1 - px < 1 - pz)
											{
												px--;
												cnx++;
											}
											else
											{
												pz--;
												cnz++;
											}
										}
									}
									else
									{
										if (pz < 1 - px)
										{
											pz++;
											cnz--;
										}
										else
										{
											px--;
											cnx++;
										}
									}
								}
							}

							tpos0 = cpos0;
							tpos1 = cpos1;
							tpos2 = cpos2;

							cpos0 = px;
							cpos1 = py;
							cpos2 = pz;

							tmin = 4;


							for (i = 0; i < 12; i++)
							{
								ttmp = (point0[i] - cpos0) * norm0[i] + (point1[i] - cpos1) * norm1[i] + (point2[i] - cpos2) * norm2[i];
								ttmp /= norm0[i] * x[0][0] + norm1[i] * x[0][1] + norm2[i] * x[0][2];

								if (ttmp > 0 && ttmp < tmin)
								{
									tmin = ttmp;
									coll = i;
								}
							}

							cnx -= norm0[coll];
							cny -= norm1[coll];
							cnz -= norm2[coll];

							cpos0 += tmin * x[0][0];
							cpos1 += tmin * x[0][1];
							cpos2 += tmin * x[0][2];

							tpos0 += tmin * x[0][0];
							tpos1 += tmin * x[0][1];
							tpos2 += tmin * x[0][2];

							cpos0 += norm0[coll];
							cpos1 += norm1[coll];
							cpos2 += norm2[coll];

							distcal += tmin;

							cface = mir[coll];

							if (cnx < 0 || cnx >= nbblocks || cny < 0 || cny >= nbblocks || cnz < 0 || cnz >= nbblocks)
							{
								currblock = 0;
							}
							else
							{
								if (cnx % 2 != 0) cnx2 = (cnx - 1) / 2;
								else cnx2 = cnx / 2;

								blockidx = cnx2 + nbblocks * cny + nbblocks * nbblocks * cnz;

								currblock = blocks[blockidx];
							}

							while (!currblock && sqrt((tpos0 - 500) * (tpos0 - 500) + (tpos1 - 500) * (tpos1 - 500) + (tpos2 - 500) * (tpos2 - 500)) < 500)
							{
								tmin = 4;


								for (i = 0; i < 12; i++)
								{
									if (i != cface) {
										ttmp = (point0[i] - cpos0) * norm0[i] + (point1[i] - cpos1) * norm1[i] + (point2[i] - cpos2) * norm2[i];
										ttmp /= norm0[i] * x[0][0] + norm1[i] * x[0][1] + norm2[i] * x[0][2];

										if (ttmp > 0 && ttmp < tmin)
										{
											tmin = ttmp;
											coll = i;
										}
									}
								}

								cnx -= norm0[coll];
								cny -= norm1[coll];
								cnz -= norm2[coll];

								cpos0 += tmin * x[0][0];
								cpos1 += tmin * x[0][1];
								cpos2 += tmin * x[0][2];

								tpos0 += tmin * x[0][0];
								tpos1 += tmin * x[0][1];
								tpos2 += tmin * x[0][2];

								distcal += tmin;

								cpos0 += norm0[coll];
								cpos1 += norm1[coll];
								cpos2 += norm2[coll];

								cface = mir[coll];

								if (cnx < 0 || cnx >= nbblocks || cny < 0 || cny >= nbblocks || cnz < 0 || cnz >= nbblocks)
								{
									currblock = false;
								}
								else
								{
									if (cnx % 2 != 0) cnx2 = (cnx - 1) / 2;
									else cnx2 = cnx / 2;
								cnx2 = cnx / 2;

								blockidx = cnx2 + nbblocks * cny + nbblocks * nbblocks * cnz;

								currblock = blocks[blockidx];
								}
							}

							if (currblock)
							{
								if (distcal>2)
								{
									buildfree = true;
								}
							}
							else
							{
								buildfree = true;
							}
						}
					}

					if (buildist < distcal) distcal = buildist;
					else distcal -= 0.1;

					cpos0 = pos0 + distcal * x[0][0];
					cpos1 = pos1 + distcal * x[0][1];
					cpos2 = pos2 + distcal * x[0][2];
					blockidx = -1;

					if (sqrt((cpos0 - 500.0) * (cpos0 - 500.0) + (cpos1 - 500.0) * (cpos1 - 500.0) + (cpos2 - 500.0) * (cpos2 - 500.0)) < 498.0)
					{
						px = fmod(cpos0, (double)1);
						py = fmod(cpos1, (double)1);
						pz = fmod(cpos2, (double)1);

						if (px < 0) px++;
						if (py < 0) py++;
						if (pz < 0) pz++;

						cnx = cpos0 - px;
						cny = cpos1 - py;
						cnz = cpos2 - pz;


						if ((cnx + cny + cnz) % 2 != 0)
						{
							if (px < py)
							{
								if (px < pz)
								{
									if (px < 1 - py)
									{
										if (px < 1 - pz)
										{
											cnx--;
										}
										else
										{
											cnz++;
										}
									}
									else
									{
										if (1 - py < 1 - pz)
										{
											cny++;
										}
										else
										{
											cnz++;
										}
									}
								}
								else
								{
									if (pz < 1 - py)
									{
										cnz--;
									}
									else
									{
										cny++;
									}
								}
							}
							else
							{
								if (py < pz)
								{
									if (py < 1 - px)
									{
										if (py < 1 - pz)
										{
											cny--;
										}
										else
										{
											cnz++;
										}
									}
									else
									{
										if (1 - px < 1 - pz)
										{
											cnx++;
										}
										else
										{
											cnz++;
										}
									}
								}
								else
								{
									if (pz < 1 - px)
									{
										cnz--;
									}
									else
									{
										cnx++;
									}
								}
							}
						}

						cnx2 = cnx / 2;

						blockidx = cnx2 + nbblocks * cny + nbblocks * nbblocks * cnz;
					}

					if (buildset)
					{
						if (!buildfree || blockidx == -1)
						{
							buildset = false;
							remidx = buildidx;
							buildidx = -1;
						}
						else if (buildidx != blockidx)
						{
							remidx = buildidx;
							buildidx = blockidx;
							addidx = buildidx;
						}
					}
					else
					{
						if (buildfree && blockidx != -1)
						{
							buildset = true;
							buildidx = blockidx;
							addidx = buildidx;
						}
					}
				}
			}
			else if (buildset == true)
			{
				buildset = false;
				remidx = buildidx;
				buildidx = -1;
			}

			if (colr == 6) col = 255;
			else col = colr + 6 * colg + 36 * colb;

			if (remidx == 498500250 ||
				remidx == 498499249 ||
				remidx == 499501250 ||
				remidx == 499500249 ||
				remidx == 499500250 ||
				remidx == 499499249 ||
				remidx == 499499250 ||
				remidx == 499498249 ||
				remidx == 500501249 ||
				remidx == 500499249 ||
				remidx == 500499250 ||
				remidx == 500500249 ||
				remidx == 500500250 ||
				remidx == 500498250 ||
				remidx == 501499250 ||
				remidx == 501500249) mmm = true;

			cudathingy(pixels, pos0, pos1, pos2, vec0, vec1, vec2, addy0, addy1, addy2, addz0, addz1, addz2,remidx,addidx,buildidx,col,nbframe);

			pixels[4 * (1280 * 360 + 640)] = 255;
			pixels[4 * (1280 * 360 + 640)+1] = 255;
			pixels[4 * (1280 * 360 + 640)+2] = 255;
			pixels[4 * (1280 * 360 + 640)+3] = 255;

			texture.update(pixels);
			sprite.setTexture(texture);
			window.draw(sprite);
			window.display();
		}
	}

	cudaExit();
	return 0;
}

void Init3(double* norm0, double* norm1, double* norm2, double* point0, double* point1, double* point2, int* mir)
{
	point0[0] = 0.5;
	point1[0] = 0.5;
	point2[0] = 1.5;

	point0[1] = 0.5;
	point1[1] = 0.5;
	point2[1] = 1.5;

	point0[2] = 0.5;
	point1[2] = 0.5;
	point2[2] = 1.5;

	point0[3] = 0.5;
	point1[3] = 0.5;
	point2[3] = 1.5;

	point0[4] = 1.5;
	point1[4] = 0.5;
	point2[4] = 0.5;

	point0[5] = 0.5;
	point1[5] = 1.5;
	point2[5] = 0.5;

	point0[6] = -0.5;
	point1[6] = 0.5;
	point2[6] = 0.5;

	point0[7] = 0.5;
	point1[7] = -0.5;
	point2[7] = 0.5;

	point0[8] = 0.5;
	point1[8] = 0.5;
	point2[8] = -0.5;

	point0[9] = 0.5;
	point1[9] = 0.5;
	point2[9] = -0.5;

	point0[10] = 0.5;
	point1[10] = 0.5;
	point2[10] = -0.5;

	point0[11] = 0.5;
	point1[11] = 0.5;
	point2[11] = -0.5;

	norm0[0] = -1;
	norm1[0] = 0;
	norm2[0] = -1;

	norm0[1] = 0;
	norm1[1] = -1;
	norm2[1] = -1;

	norm0[2] = 1;
	norm1[2] = 0;
	norm2[2] = -1;

	norm0[3] = 0;
	norm1[3] = 1;
	norm2[3] = -1;

	norm0[4] = -1;
	norm1[4] = -1;
	norm2[4] = 0;

	norm0[5] = 1;
	norm1[5] = -1;
	norm2[5] = 0;

	norm0[6] = 1;
	norm1[6] = 1;
	norm2[6] = 0;

	norm0[7] = -1;
	norm1[7] = 1;
	norm2[7] = 0;

	norm0[8] = -1;
	norm1[8] = 0;
	norm2[8] = 1;

	norm0[9] = 0;
	norm1[9] = -1;
	norm2[9] = 1;

	norm0[10] = 1;
	norm1[10] = 0;
	norm2[10] = 1;

	norm0[11] = 0;
	norm1[11] = 1;
	norm2[11] = 1;

	mir[0] = 10;
	mir[1] = 11;
	mir[2] = 8;
	mir[3] = 9;
	mir[4] = 6;
	mir[5] = 7;
	mir[6] = 4;
	mir[7] = 5;
	mir[8] = 2;
	mir[9] = 3;
	mir[10] = 0;
	mir[11] = 1;
}
